import os
import numpy as np
import collections
import argparse
import common.libbgmg

def check_input_file(args, argname, chri=None):
    arg_dict = vars(args)
    if (argname in arg_dict) and arg_dict[argname]:
        argval = arg_dict[argname]
        if chri: argval=argval.replace('@', str(chri))
        if not os.path.isfile(argval): raise ValueError("Input file --{a} does not exist: {f}".format(a=argname, f=argval))

def fix_chr2use_arg(args, argname="chr2use"):
    arg_dict = vars(args)
    if argname not in arg_dict : return
    if not arg_dict[argname]: arg_dict[argname] = []; return
    chr2use_arg = arg_dict[argname]
    chr2use = []
    if chr2use_arg is not None:
        for a in chr2use_arg.split(","):
            if "-" in a:
                start, end = [int(x) for x in a.split("-")]
                chr2use += [str(x) for x in range(start, end+1)]
            else:
                chr2use.append(a.strip())
        if np.any([not x.isdigit() for x in chr2use]): raise ValueError('Chromosome labels must be integer')
    arg_dict[argname] = [int(x) for x in chr2use]

def make_ranges(args_exclude_ranges):
    # Interpret --exclude-ranges input
    ChromosomeRange = collections.namedtuple('ChromosomeRange', ['chr', 'from_bp', 'to_bp'])
    exclude_ranges = []
    if args_exclude_ranges is not None:
        for exclude_range in args_exclude_ranges:
            try:
                exclude_range=exclude_range.strip().upper().replace('CHR','')
                if exclude_range == '[]': continue
                if exclude_range == 'MHC': exclude_range='6:25-35MB'
                if exclude_range == 'APOE': exclude_range = '19:45-46MB'
                if exclude_range.endswith('KB'): factor=1e3
                elif exclude_range.endswith('MB'): factor=1e6
                else: factor=1
                exclude_range=exclude_range.replace('KB', '').replace('MB', '')
                chr_from_to = [int(x) for x in exclude_range.replace(':', ' ').replace('-', ' ').split()[:3]]
                chr_from_to[1] *= factor
                chr_from_to[2] *= factor
                range = ChromosomeRange._make(chr_from_to)
            except Exception as e:
                raise(ValueError('Unable to interpret exclude range "{}", chr:from-to format is expected.'.format(exclude_range)))
            exclude_ranges.append(range)
    return exclude_ranges

def fix_and_validate_args(args):
    # here is what you need to know about how --seed  is handled:
    # 1. If it optional. When not specified, we will use np.rando.randint to generate a random seed (the the line below this comment)
    # 2. There are several places in this code where we generate random numbers using numpy library. Therefore, we use args.seed to initialize np.random.seed.
    # 3. We pass args.seed to MiXeR using set_options("seed", args.seed), see convert_args_to_libbgmg_options() function.
    #    This will ensure that random prunning uses consistent seed.
    # 4. In ADAM algorithm, each cost function evaluation starts with it's own seed (to respect the stochastic nature of the algorithm)
    #    This is done in line "libbgmg.set_option('seed', np.random.randint(np.iinfo(np.int32).max))".
    # 5. --diffevo-fast-repeats also rely on seed to ensure different path of the differential evolution algorihtm.
    if 'seed' in args:
        if args.seed is None: args.seed = np.random.randint(np.iinfo(np.int32).max)
        np.random.seed(args.seed)

    fix_chr2use_arg(args, argname="chr2use")
    if 'trait1_file' in args:
        if '@' in args.trait1_file:
            for chri in args.chr2use:
                check_input_file(args, 'trait1_file', chri)
        else:
            check_input_file(args, 'trait1_file')
    check_input_file(args, 'trait2_file')
    check_input_file(args, 'trait1_params_file')
    check_input_file(args, 'trait2_params_file')
    check_input_file(args, 'load-params-file')

    if ('analysis' in args) and (args.analysis == 'fit2'):
        if not args.load_params_file:
            if not args.trait1_params_file: raise ValueError('--trait1-params-file or --load-params-file is required ')
            if not args.trait2_params_file: raise ValueError('--trait2-params-file or --load-params-file is required ')

    arg_dict = vars(args)
    for chri in arg_dict["chr2use"]:
        check_input_file(args, 'bim-file', chri)
        check_input_file(args, 'ld-file', chri)
        check_input_file(args, 'annot-file', chri)

# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string=None):
        with values as f:
            contents = '\n'.join([x for x in f.readlines() if (not x.strip().startswith('#'))])

        data = parser.parse_args(contents.split(), namespace=namespace)
        for k, v in vars(data).items():
            if v and k != option_string.lstrip('-'):
                setattr(namespace, k, v)

def initialize_libbgmg_logging(args, header):
    libbgmg = common.libbgmg.LibBgmg(args.lib)
    log_fname = args.log if args.log else args.out + '.log'
    if os.path.exists(log_fname): os.system(f'rm {log_fname}')
    libbgmg.init_log(log_fname)
    libbgmg.log_message(header)
    libbgmg.dispose()
    return libbgmg

def convert_args_to_libbgmg_options(args):
    libbgmg_options = {
        'r2min': args.r2min if ('r2min' in args) else None,
        'kmax': args.kmax[0] if ('kmax' in args) else None, 
        'kmax_pdf': args.kmax_pdf[0] if ('kmax_pdf' in args) else None,
        'threads': args.threads[0] if ('threads' in args) else None,
        'seed': args.seed if ('seed' in args) else None,
        'cubature_rel_error': args.cubature_rel_error if ('cubature_rel_error' in args) else None,
        'cubature_max_evals': args.cubature_max_evals if ('cubature_max_evals' in args) else None,
        'z1max': args.z1max if ('z1max' in args) else None,
        'z2max': args.z2max if ('z2max' in args) else None, 
    }
    return [(k, v) for k, v in libbgmg_options.items() if v is not None ]
