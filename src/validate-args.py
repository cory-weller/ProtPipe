#!/usr/bin/env python3
if __name__ == '__main__':

    #### VERSION ##################################################################################
    singularity_image = "src/diann-1.8.1.sif"
    singularity_image_md5 = '35644c1d7217f0c65727b8fb9c8bfaae'
    diann_version = "1.8.1"
    singularity_version = "3.x"
    git_repo = "https://www.github.com/cory-weller/ProtPipe"
    authors = ['Cory Weller']

    #### LOGGING ###################################################################################
    import logging
    class ExitOnExceptionHandler(logging.StreamHandler):
        '''logger class to exit on error or critical'''
        def emit(self, record):
            super().emit(record)
            if record.levelno in (logging.ERROR, logging.CRITICAL):
                raise SystemExit(-1)

    logger = logging.getLogger(None)
    logging.basicConfig(handlers=[ExitOnExceptionHandler()],
                        level=logging.DEBUG,
                        format='%(message)s'
    )
    logger.arg_errors = []
    logger.line_errors = []
    logger.permission_errors = []


    #### MODULES ####################################################################################
    import subprocess
    import os
    import sys
    import re
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--dry',
                        default=False,
                        action="store_true",
                        help='Validate arguments and file permissions, but do not run DIA-NN')
    parser.add_argument('--cfg', 
                        default='config.txt',
                        nargs='?',
                        type=str,
                        metavar='<config.txt>',
                        help='File name to import DIA-NN configuration. Default: config.txt')

    args = parser.parse_args()

    # if len(sys.argv) < 4:
    #     logger.error("ERROR: Provide config file and singularity image file")
    
    # if len(sys.argv) > 4:
    #     logger.error("ERROR: too many arguments given. Only specify one config file.")
    
    # if len(sys.argv) == 4:
    #     config_filename = sys.argv[1]
    #     logger.info(f"INFO: Using config file {config_filename}")
    #     diann_sif_filename = sys.argv[2]
    #     logger.info(f"INFO: Using DIA-NN singularity image file {diann_sif_filename}")
    #     dryrun = sys.argv[3]
    #     if dryrun:
    #         logger.info("INFO: --dry-run indicated, halting after checks, but before running DIA-NN")


    #### FUNCTIONS #################################################################################
    def collect_arguments(arg, val, lineNo, arg_collector, option_collector):
        arg = '--'+arg
        if not (arg in valid_args or arg in valid_options):
            logger.arg_errors.append(f"  Line {lineNo}: {arg} is not a recognized option")
        if arg in valid_args:
            if val is None:
                logger.arg_errors.append(f"  Line {lineNo}: {arg} requires a value along with it")
            try:
                arg_collector[arg].append(val)
                logger.info(f"INFO: {arg} already exists, adding additional value {val}")
            except KeyError:
                arg_collector[arg]=[val]
        if arg in valid_options:
            if val is not None:
                logger.arg_errors.append(f"  Line {lineNo}: {arg} does not take arguments, ignoring {val}")
            elif arg not in option_collector:
                option_collector.append(arg)
    
    def basedir(path):
        path = re.split(string=path, pattern='/+')
        path = '/'.join(path[:-1])
        return path

    def read_config(filename, logger):
        try:
            with open(filename, 'r') as infile:
                for line in infile:
                    yield line
        except FileNotFoundError:
            logger.arg_errors.append(f"ERROR: config file {filename} not found")

    #### DEFINITIONS ###############################################################################
    collected_args = {}
    collected_options = []

    valid_args = [
    '--f',
    '--out',
    '--fasta',
    '--dir',
    '--threads',
    '--out-lib',
    '--qvalue',
    '--min-fr-mz',
    '--max-fr-mz',
    '--prefix',
    '--decoy-channel',
    '--ext',
    '--fixed-mod',
    '--im-window',
    '--im-window-factor',
    '--learn-lib',
    '--lib',
    '--channels',
    '--cut',
    '--lib-fixed-mod',
    '--max-pep-len',
    '--max-pr-charge',
    '--min-pep-len',
    '--min-pr-charge',
    '--min-pr-mz',
    '--max-pr-mz',
    '--no-cut-after-mod',
    '--library-headers',
    '--matrix-ch-qvalue',
    '--matrix-qvalue',
    '--matrix-tr-qvalue',
    '--sptxt-acc',
    '--mass-acc',
    '--mass-acc-cal',
    '--mass-acc-ms1',
    '--max-fr',
    '--missed-cleavages',
    '--min-fr',
    '--min-peak',
    '--mod',
    '--monitor-mod',
    '--pg-level',
    '--quant-fr',
    '--target-fr',
    '--verbose',
    '--var-mod',
    '--var-mods',
    '--vis',
    '--window',]

    required_arg_defaults = {
        '--threads' : 4,
    }

    valid_options=[
    '--met-excision',
    '--matrix-spec-q',
    '--mbr-fix-settings',
    '--peak-center',
    '--peak-height',
    '--peak-translation',
    '--nn-single-seq',
    '--mod-only',
    '--predictor',
    '--reanalyse',
    '--reannotate',
    '--regular-swath',
    '--relaxed-prot-inf',
    '--report-lib-info',
    '--restrict-fr',
    '--scanning-swath',
    '--smart-profiling',
    '--species-genes',
    '--strip-unknown-mods',
    '--no-lib-filter',
    '--tims-skip-errors',
    '--use-quant',
    '--no-calibration',
    '--no-decoy-channel',
    '--no-fr-selection',
    '--no-ifs-removal',
    '--no-im-window',
    '--no-isotopes',
    '--no-main-report',
    '--no-maxlfq',
    '--no-norm',
    '--no-prot-inf',
    '--no-quant-files',
    '--no-rt-window',
    '--no-stats',
    '--no-stratification',
    '--no-swissprot',
    '--original-mods',
    '--duplicate-proteins',
    '--il-eq',
    '--quick-mass-acc',
    '--semi',
    '--convert',
    '--gen-spec-lib',
    '--matrices',
    '--out-lib-copy',
    '--out-measured-rt',
    '--clear-mods',
    '--compact-report',
    '--dl-no-im',
    '--dl-no-rt',
    '--exact-fdr',
    '--force-swissprot',
    '--full-unimod',
    '--gen-fr-restriction',
    '--global-mass-cal',
    '--global-norm',
    '--individual-mass-acc',
    '--individual-reports',
    '--individual-windows',
    '--int-removal',
    '--fasta-search',]

    required_option_defaults = {

    }

    #### RUN #######################################################################################

    # Iterate over config_filename and collect args, logging if malformed or invalid
    # Lines commented out will be ignored
    # if 'IGNORE' is encountered, break out of loop
    DIA_NN_args = read_config(args.cfg, logger)
    for i,rawline in enumerate(DIA_NN_args):
        if rawline.startswith('IGNORE'):
            break
        line = rawline.split('#')[0].strip('\n .-').split()
        if len(line)==0:        # ignore blank or commented out lines
            continue
        if len(line) > 2:       # possibly malformed
            logger.line_errors.append(f"  Line {i+1}: {rawline.rstrip()}")
        argname=line[0]
        if len(line) == 2:
            argval=line[1]
        else:
            argval = None
        collect_arguments(argname, argval, i+1, collected_args, collected_options)

    # Print any errors after parsing whole config file
    if len(logger.line_errors) > 0:
        print("ERROR: Check malformed lines and try again:")
        print('\n'.join(logger.line_errors)+'\n')

    if len(logger.arg_errors) > 0:
        print("ERROR: Check argument(s):")
        print('\n'.join(logger.arg_errors)+'\n')


    # Check output directory is writable
    if '--out' in collected_args:
        outdirs = collected_args['--out']   # entire list
        outdir = outdirs[-1]                # last item in list, regardless of length
        if len(outdirs) > 1:
            logger.warning(f"WARNING: multiple --out arguments given. Only using {outdir}")
            collected_args['--out'] = [outdir]
        if not any([outdir.endswith(i) for i in ['txt','tsv']]):
            msg = f"ERROR: --out must specify an output file ending in .txt or .tsv\n"
            msg += f"You specified: --out {outdir}"
            print(msg)
            logger.permission_errors.append(msg)
        outdir = basedir(outdir)
        if not os.path.isdir(outdir):
            logger.info(f"INFO: trying to create --out directory {outdir}")
            if not args.dry:
                try: os.makedirs(outdir)
                except PermissionError:
                    msg = f"ERROR: permission denied to create directory {outdir}"
                    logger.permission_errors.append(msg)



    total_errors = len(logger.arg_errors + logger.line_errors + logger.permission_errors)
    if total_errors > 0:
        logger.error(f"Exiting due to {total_errors} problem(s)")
    elif total_errors == 0:
        logger.info("INFO: No troubling errors found")
        # diann_args = format_args(collected_args, collected_options)
        #print(collected_args)
        final_args = ['singularity', 'exec', '--cleanenv', '-H', os.getcwd(), singularity_image, 'diann']
        for key in collected_args:
            for value in collected_args[key]:
                final_args.append(f"{key} {value}")
        for value in collected_options:
            final_args.append(value)
        logger.info('INFO: command to be called:\n'+' '.join(final_args))
        if args.dry:
            logger.info("INFO: Stopping before running DIA-NN because --dry-run was specified")
        else:
            subprocess.run(final_args)