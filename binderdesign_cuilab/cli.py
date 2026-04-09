import argparse
import sys


def _project_main(argv):
    from . import project

    parser = project.build_parser()
    args = parser.parse_args(argv)
    args.func(args)


def _table_main(argv):
    from . import table

    parser = table.build_parser()
    args = parser.parse_args(argv)
    args.func(args)


def task_main(argv):
    parser = argparse.ArgumentParser(prog="binderdesign-cuilab task")
    subparsers = parser.add_subparsers(dest="task", required=True)

    rfd = subparsers.add_parser("rfdiffusion")
    rfd.add_argument("--target-pdb", required=True)
    rfd.add_argument("--config-file", required=True)
    rfd.add_argument("--output-root", required=True)
    rfd.add_argument("--num-designs", type=int, required=True)
    rfd.add_argument("--index", type=int, required=True)

    qc = subparsers.add_parser("light-qc")
    qc.add_argument("--manifest", required=True)
    qc.add_argument("--index", type=int, required=True)
    qc.add_argument("--template-pdb", required=True)
    qc.add_argument("--output-root", required=True)
    qc.add_argument("--result-dir", required=True)

    geom = subparsers.add_parser("geometry")
    geom.add_argument("--manifest", required=True)
    geom.add_argument("--index", type=int, required=True)
    geom.add_argument("--result-dir", required=True)
    geom.add_argument("--hotspots", required=True)
    geom.add_argument("--binder-length-min", type=int, required=True)
    geom.add_argument("--binder-length-max", type=int, required=True)
    geom.add_argument("--pdb-column", default="restored_pdb")

    mpnn = subparsers.add_parser("proteinmpnn")
    mpnn.add_argument("--manifest", required=True)
    mpnn.add_argument("--index", type=int, required=True)
    mpnn.add_argument("--output-root", required=True)
    mpnn.add_argument("--pdb-column", default="restored_pdb")

    af = subparsers.add_parser("alphafold", aliases=["af2", "af3"])
    af.add_argument("--manifest", required=True)
    af.add_argument("--index", type=int, required=True)
    af.add_argument("--output-root", required=True)
    af.add_argument("--backend", choices=["af2", "af3"])

    ia = subparsers.add_parser("rosetta-interface")
    ia.add_argument("--manifest", required=True)
    ia.add_argument("--index", type=int, required=True)
    ia.add_argument("--work-root", required=True)
    ia.add_argument("--result-dir", required=True)
    ia.add_argument("--pdb-column", default="restored_pdb")

    fr = subparsers.add_parser("rosetta-fastrelax")
    fr.add_argument("--manifest", required=True)
    fr.add_argument("--index", type=int, required=True)
    fr.add_argument("--work-root", required=True)
    fr.add_argument("--raw-dir", required=True)
    fr.add_argument("--pdb-column", default="restored_pdb")

    ddg = subparsers.add_parser("rosetta-ddg")
    ddg.add_argument("--manifest", required=True)
    ddg.add_argument("--index", type=int, required=True)
    ddg.add_argument("--work-root", required=True)
    ddg.add_argument("--result-dir", required=True)
    ddg.add_argument("--pdb-column", default="restored_pdb")

    agg = subparsers.add_parser("aggregate-robustness")
    agg.add_argument("--raw-dir", required=True)
    agg.add_argument("--output", required=True)

    args = parser.parse_args(argv)

    if args.task == "rfdiffusion":
        from . import tasks

        tasks.run_rfdiffusion_task(args.target_pdb, args.config_file, args.output_root, args.num_designs, args.index)
        return
    if args.task == "light-qc":
        from . import tasks

        tasks.run_light_qc_task(args.manifest, args.index, args.template_pdb, args.output_root, args.result_dir)
        return
    if args.task == "geometry":
        from . import tasks

        tasks.run_geometry_task(
            args.manifest,
            args.index,
            args.result_dir,
            args.hotspots,
            args.binder_length_min,
            args.binder_length_max,
            args.pdb_column,
        )
        return
    if args.task == "proteinmpnn":
        from . import tasks

        tasks.run_proteinmpnn_task(args.manifest, args.index, args.output_root, args.pdb_column)
        return
    if args.task in {"alphafold", "af2", "af3"}:
        from . import tasks

        backend = args.backend or ("af3" if args.task == "af3" else "af2" if args.task == "af2" else None)
        tasks.run_alphafold_task(args.manifest, args.index, args.output_root, backend=backend)
        return
    if args.task == "rosetta-interface":
        from . import tasks

        tasks.run_rosetta_interface_task(args.manifest, args.index, args.work_root, args.result_dir, args.pdb_column)
        return
    if args.task == "rosetta-fastrelax":
        from . import tasks

        tasks.run_rosetta_fastrelax_task(args.manifest, args.index, args.work_root, args.raw_dir, args.pdb_column)
        return
    if args.task == "rosetta-ddg":
        from . import tasks

        tasks.run_rosetta_ddg_task(args.manifest, args.index, args.work_root, args.result_dir, args.pdb_column)
        return
    if args.task == "aggregate-robustness":
        from . import tasks

        tasks.aggregate_robustness(args.raw_dir, args.output)
        return


def build_parser():
    parser = argparse.ArgumentParser(prog="binderdesign-cuilab", description="BinderDesign-CUILab end-to-end binder design workflow library")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init-project", help="Initialize a BinderDesign-CUILab project")
    init.add_argument("remainder", nargs=argparse.REMAINDER)

    submit = subparsers.add_parser("submit-marlowe", help="Submit a BinderDesign-CUILab project on Marlowe")
    submit.add_argument("remainder", nargs=argparse.REMAINDER)

    run = subparsers.add_parser("run-project", help="Initialize and optionally submit a project")
    run.add_argument("remainder", nargs=argparse.REMAINDER)

    project_group = subparsers.add_parser("project", help="Project orchestration subcommands")
    project_group.add_argument("remainder", nargs=argparse.REMAINDER)

    table_group = subparsers.add_parser("table", help="Master-table operations")
    table_group.add_argument("remainder", nargs=argparse.REMAINDER)

    task_group = subparsers.add_parser("task", help="Single-stage workers used by Slurm")
    task_group.add_argument("remainder", nargs=argparse.REMAINDER)

    return parser


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv:
        if argv[0] == "project":
            _project_main(argv[1:])
            return
        if argv[0] == "table":
            _table_main(argv[1:])
            return
        if argv[0] == "task":
            task_main(argv[1:])
            return
        if argv[0] in {"init-project", "run-project", "submit-marlowe"}:
            _project_main(argv)
            return
    parser = build_parser()
    if not argv:
        parser.print_help()
        return
    parser.parse_args(argv)
    parser.print_help()


if __name__ == "__main__":
    main()
