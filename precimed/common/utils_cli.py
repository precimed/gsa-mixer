import os
import numpy as np
import collections
import argparse
import common.libbgmg

GO_EXTEND_BP_DEFAULT = 10000

def parent_parser_add_argument_out_lib_log(parent_parser):
    parent_parser.add_argument('--argsfile', type=open, action=common.utils_cli.LoadFromFile, default=None, help="file with additional command-line arguments, e.g. those that are across all of your mixer.py runs (--lib, --bim-file and --ld-file)")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin")
    parent_parser.add_argument("--log", type=str, default=None, help="file to output log, defaults to <out>.log")

def parser_add_argument_bim_file(parser):
    parser.add_argument("--bim-file", type=str, default=None, help="plink bim file (required argument); "
        "defines the reference set of SNPs used for the analysis. "
        "Marker names must not have duplicated entries. "
        "May contain symbol '@', which will be replaced by an actual chromosome label. ")

def parser_add_argument_ld_file(parser):
    parser.add_argument("--ld-file", type=str, default=None, help="file with linkage disequilibrium information, "
        "generated via 'mixer.py ld' command (required argument); "
        "may contain symbol '@', similarly to --bim-file argument. ")

def parser_add_argument_chr2use(parser):
    parser.add_argument("--chr2use", type=str, default="1-22", help="chromosome labels to use (default: %(default)s); "
        "chromosome must be labeled by integer, i.e. X and Y are not acceptable; example of valid arguments: '1,2,3' or '1-4,12,16-20'")

def parser_add_argument_extract_and_exclude(parser):
    parser.add_argument('--extract', type=str, default="", help="(optional) File with variants to include in the analysis. By default, all variants are included. This applies to GWAS tag SNPs, however LD is still computed towards the entire reference provided by --bim-file. This applies before --exclude, so a variant listed in --exclude won't be re-introduced if it's also present in --extract list.")
    parser.add_argument('--exclude', type=str, default="", help="(optional) File with variants to exclude from the analysis.")
    parser.add_argument('--exclude-ranges', type=str, default=['MHC'], nargs='+',
        help='(default: %(default)s) exclude SNPs in ranges of base pair position; '
        'the syntax is chr:from-to, for example 6:25000000-35000000; '
        'multiple regions can be excluded; '
        '"chr" prefix prior to chromosome label, as well as KB and MB suffices are allowed, e.g. chr6:25-35MB is a valid exclusion range. '
        'Some special case regions are also supported, for example "--exclude-ranges MHC APOE".'
        'To overwrite the default, pass "--exclude-ranges []".')

def parser_add_argument_use_complete_tag_indices(parser):
    parser.add_argument('--use-complete-tag-indices', default=False, action="store_true", help=argparse.SUPPRESS) #"advanced option (expert use only); MiXeR software has a notion of 'tag SNPs', which are the subset of the SNPs in --bim-file with well defined GWAS z-score(s); this saves memory used by LD r2 storage, which became a rectangular matrix (LD r2 between two SNPs from the --bim-file is stored only if one of the SNPs is a tag SNP). Setting this flag enforces all --bim-file SNPs to be treated as tag SNPs, some will have undefined GWAS z-score(s), and the complete LD matrix will be save in memory.")

def parser_add_argument_allow_albiguous_snps(parser):
    parser.add_argument('--allow-ambiguous-snps', default=False, action="store_true", help="advanced option (expert use only); a flag allowing to include A/T and C/G SNPs in fit procedure.")

def parser_add_argument_savelib_and_loadlib(parser):
    parser.add_argument("--savelib-file", type=str, default=None, help="optional path to save the state of the libbgmg object (for faster initialization with --loadlib-file). It is highly recommended to use --use-complete-tag-indices together with --savelib-file.")
    parser.add_argument("--loadlib-file", type=str, default=None, help="optional path to load libbgmg object; "
        "there are several limitations if you use this option: "
        "--bim-file argument must be the same as what was used with --savelib-file; "
        "for plsa analysis --chr2use can be different as each chromosome has its own .bin file; "
        "--ld-file and --use-complete-tag-indices arguments will be ignored.")

def parser_add_argument_annot_file(parser, annot_file_test=False):
    parser.add_argument("--annot-file", type=str, default=None, help="(optional) path to binary annotations in LD score regression format, i.e. <path>/baseline.@.annot.gz for fitting enrichment model model. This must include the first column with all ones ('base' annotation category covering the entire genome).")
    if annot_file_test:
        parser.add_argument("--annot-file-test", type=str, default=None, help="(optional) path to binary annotations in LD score regression format, i.e. <path>/baseline.@.annot.gz for evaluating enrichment.  If provided, enrichment scores computed for --annot-file-test will be saved to <out>.annot_test_enrich.csv.")

def parser_add_argument_go_file(parser, go_file_test=False):
    # --go-file and --go-file-test must contain columns 'GO', 'GENE', 'CHR', 'FROM', 'TO'
    parser.add_argument("--go-file", type=str, default=None, help="(optional) path to GO antology file for fitting enrichment model model. The format is described in the documentation. 'base' category that covers entire genome will be added automatically.")
    parser.add_argument("--go-extend-bp", type=int, default=GO_EXTEND_BP_DEFAULT, help="extends each gene by this many base pairs, defining a symmetric window up and downstream (default: %(default)s)")
    if go_file_test:
        parser.add_argument("--go-file-test", type=str, default=None, help="(optional) path to an additional  GO antology file for evaluating enrichment, same convention as --go-file. If provided, enrichment scores computed for --go-test-file will be saved to <out>.go_test_enrich.csv")

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
