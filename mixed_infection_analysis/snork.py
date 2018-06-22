#!/usr/bin/env python

import argparse
import logging
import os
import random
import sys
import time

import snorkversion

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

_errorargsinvalid = 3

def Assume_Bin_Is_Same_Dir_As_Executable():
    bin1 = '/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-2] + ['bin'])
    bin2 = os.path.dirname(os.path.realpath(sys.argv[0]))
    sys.path.insert(0, bin1)
    sys.path.insert(0, bin2)
    _bin = bin1
    return _bin

def Parse_Bin_From_Arguments():
    _bin = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('-bin'):
                if '=' in sys.argv[i]:
                    _bin = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                    break
                else:
                    _bin = os.path.realpath(os.path.expandvars(sys.argv[i+1]))
                    break
    except:
        _bin = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if _bin is None:
        _bin = os.path.dirname(os.path.realpath(os.path.expandvars(sys.argv[0])))
    if _bin and os.path.exists(_bin):
        try:
            sys.path.insert(0, _bin)
        except:
            sys.stderr.write('Error: Failed to import snorklib module\n')
            sys.exit(_errorargsinvalid)
    else:
        _bin = Assume_Bin_Is_Same_Dir_As_Executable()
    return _bin

def Parse_Profile_From_Arguments(bindir):
    _profile = None
    try:
        for i in range(1, len(sys.argv)):
            if sys.argv[i].startswith('--profile'):
                if '=' in sys.argv[i]:
                    _profile = os.path.realpath(os.path.expandvars('='.join(sys.argv[i].split('=')[1:])))
                else:
                    _profile = os.path.realpath(os.path.expandvars(sys.argv[i+1]))
    except:
        if os.path.exists(bindir):
            _profile = os.path.join(bindir, 'snork.profile')
        else:
            _profile = 'snork.profile'
    if _profile is None:
        _profile = os.path.join(bindir, 'snork.profile')
    if os.path.exists(_profile):
        with open (_profile, 'r') as in_fp:
            for line in in_fp:
                if not line.startswith('export'):
                    continue
                info = line.strip().split(' ')[1].split('=')
                if len(info) == 2:
                    var, val = info
                    sys.stderr.write('Setting {0}={1}\n'.format(var, val))
                    os.environ[var] = val
    return _profile

def Parse_Bin_And_Profile():
    _bin = None
    _profile = None
    if len(sys.argv) == 1:
        _bin = Assume_Bin_Is_Same_Dir_As_Executable()
        import snorklib
    else:
        _bin = Parse_Bin_From_Arguments()
        _profile = Parse_Profile_From_Arguments(_bin)
    return _bin, _profile

def run_subtool(parser, args, P, mylogger, myhandler):
    if args.command == 'splitpops':
        import splitpops as submodule
    else:
        return
    submodule.main(parser, args, P, mylogger, myhandler, sys.argv)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)

def Create_Program_Parser_And_Subparsers():
    global parser
    global subparsers
    parser = argparse.ArgumentParser(prog=snorkversion.__program__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '-version', help='Print version', action='version',
        version=snorkversion.__program__+' version '+snorkversion.__version__)
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command',
        parser_class=ArgumentParserWithDefaults)

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def Create_Splitpops_Subparser():
    global subparsers
    p = subparsers.add_parser('splitpops', help='Split target organism reads into populations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-bin', dest='bin', metavar='DIR', required=False,
        default='.', help='Scripts dir')
    p.add_argument('-profile', dest='profile', metavar='FILE', required=False,
        default=None, help='Shell statements to set environment variables')
    p.add_argument('-config', dest='config', metavar='FILE', required=False,
        default=None, help='Run file of config options')
    p.add_argument('-orgid', dest='orgid', metavar='STR', required=True,
        default='UNKNOWN', help='Organism identifier')
    p.add_argument('-dataid', dest='dataid', metavar='STR', required=True,
        default='UNKNOWN', help='Primary identifier for outfiles')
    p.add_argument('-samplename', dest='samplename', metavar='STR', required=True,
        default='UNKNOWN', help='Sample name')
    p.add_argument('-bampath', dest='bampath', metavar='FILE', required=False,
        default=None, help='Paired-end reads from sample')
    p.add_argument('-targetrefid', dest='targetrefid', metavar='STR', required=True,
        default='UNKNOWN', help='ID for target references of known genotype class')
    p.add_argument('-targetrefpath', dest='targetrefpath', metavar='FILE', required=True,
        default='UNKNOWN', help='Path of target refs each with contigid of format SUBGT_XXXX ' \
        ' and seq being the genome of each strain (FASTA)')
    p.add_argument('-mintargetpoppct', dest='mintargetpoppct', metavar='FILE', required=False, type=float,
        default=1.0, help='Do not report target populations with <threshold%% reads')
    p.add_argument('-outdir', dest='outdir', metavar='DIR', required=True,
        default='.', help='Output directory')
    p.add_argument('-logdir', dest='logdir', metavar='DIR', required=True,
        default='.', help='Logging directory')
    p.add_argument('-dryrun', dest='dryrun', metavar='BOOL', required=False, type=str2bool,
        default=False, help='Print tasks without running')
    p.add_argument('-overwrite', dest='overwrite', metavar='BOOL', required=False, type=str2bool,
        default=False, help='Overwrite outfiles')
    p.add_argument('-deleteints', dest='deleteints', metavar='BOOL', required=False, type=str2bool,
        default=True, help='Delete intermediate outfiles')
    p.add_argument('-verbosity', dest='verbosity', metavar='STR', required=False,
        default='INFO', help='Diagnostics level [DEBUG|INFO|ERROR]')
    p.set_defaults(func=run_subtool)
    return p

def Setup_Diagnostic_Logging(args, _P):
    now = time.localtime()
    logdir = None
    try:
         logdir = args.logdir
    except:
        sys.stderr.write('Error: Failed to set log dir\n')
        sys.exit(_P.err_code('ErrorDirCreate'))
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    logfile = '{0}-snork-{1:03d}.log'.format(
        time.strftime('%Y%m%d-%H%M%S', now), random.randint(1, 100))
    logpath = os.path.join(logdir, logfile)
    handler = logging.FileHandler(logpath)
    handler.setLevel(logging.INFO)
    loggingverbosity = args.verbosity.upper()
    if loggingverbosity == 'ERROR':
        logger.setLevel(logging.ERROR)
        handler.setLevel(logging.ERROR)
    elif loggingverbosity == 'DEBUG':
        logger.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(levelname)s : %(message)s', '%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info('Program: '+snorkversion.__program__+' version '+snorkversion.__version__+' ('+snorkversion.__description__+')')
    logger.info('Command: '+' '.join(sys.argv))
    return handler, logger

def main():
  # Parse bin dir and profile path to get defaults from config file
    _bin, _profile = Parse_Bin_And_Profile()
    import snorklib
    _P = snorklib.snorklib()
  # Parse command-line
    Create_Program_Parser_And_Subparsers()
    p01 = Create_Splitpops_Subparser()
    args = parser.parse_args()
  # Fix type of args
  # Overwrite the command-line values for bin and profile
    if not args.bin or args.bin is None:
        args.bin = _bin
    if not args.profile or args.profile is None:
        args.profile = _profile
  # Read the -config file
    if not _P.config_read(args.config):
        sys.stderr.write('Error: Failed to read conf file ({0})\n'.format(args.config))
        sys.exit(_P.err_code('ErrorReadingData'))
  # Set logging verbosity and output file
    handler, logger = Setup_Diagnostic_Logging(args, _P)
  # Call the selected sub-command, trap signals but ignore SIGPIPEs (32)
    try:
        args.func(parser, args, _P, logger, handler)
    except IOError, e :
        if e.errno != 32:
            logger.info('Received kill signal')
            raise

if __name__ == '__main__':
    main()
