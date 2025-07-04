"""
Command line interface for Palace Python package.
"""

import argparse
import os
import sys

from .core import PalaceSolver
from .utils import create_basic_config, get_example_configs, load_config, save_config


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Palace Python Interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  palace-py run config.json                 # Run simulation
  palace-py run config.json -np 4           # Run with 4 MPI processes
  palace-py run config.json -o results/     # Specify output directory
  palace-py create-config mesh.msh          # Create basic config
  palace-py list-examples                   # List available examples
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Run command
    run_parser = subparsers.add_parser("run", help="Run Palace simulation")
    run_parser.add_argument("config", help="Configuration file")
    run_parser.add_argument("-o", "--output-dir", help="Output directory")
    run_parser.add_argument(
        "-np", "--num-procs", type=int, help="Number of MPI processes"
    )
    run_parser.add_argument("--palace-exe", help="Path to palace executable")
    run_parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output"
    )

    # Create config command
    config_parser = subparsers.add_parser(
        "create-config", help="Create basic configuration"
    )
    config_parser.add_argument("mesh_file", help="Mesh file path")
    config_parser.add_argument(
        "-t",
        "--type",
        choices=["Eigenmode", "Driven", "Transient"],
        default="Eigenmode",
        help="Problem type",
    )
    config_parser.add_argument("-o", "--output", help="Output config file name")

    # List examples command
    subparsers.add_parser("list-examples", help="List available example configurations")

    # Validate command
    validate_parser = subparsers.add_parser(
        "validate", help="Validate configuration file"
    )
    validate_parser.add_argument("config", help="Configuration file to validate")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    try:
        if args.command == "run":
            return run_simulation(args)
        elif args.command == "create-config":
            return create_config_file(args)
        elif args.command == "list-examples":
            return list_examples()
        elif args.command == "validate":
            return validate_config(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


def run_simulation(args):
    """Run Palace simulation."""
    solver = PalaceSolver(args.palace_exe)

    print(f"Running Palace simulation: {args.config}")
    if args.num_procs:
        print(f"  MPI processes: {args.num_procs}")
    if args.output_dir:
        print(f"  Output directory: {args.output_dir}")

    result = solver.run(
        args.config, output_dir=args.output_dir, num_procs=args.num_procs
    )

    if args.verbose or result.returncode != 0:
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)

    if result.returncode == 0:
        print("Simulation completed successfully!")
    else:
        print(f"Simulation failed with exit code {result.returncode}")

    return result.returncode


def create_config_file(args):
    """Create basic configuration file."""
    if not os.path.exists(args.mesh_file):
        print(f"Error: Mesh file not found: {args.mesh_file}", file=sys.stderr)
        return 1

    config = create_basic_config(args.mesh_file, args.type)

    output_file = args.output or f"{args.type.lower()}_config.json"
    save_config(config, output_file)

    print(f"Created configuration file: {output_file}")
    print(f"  Problem type: {args.type}")
    print(f"  Mesh file: {args.mesh_file}")

    return 0


def list_examples():
    """List available example configurations."""
    examples = get_example_configs()

    if not examples:
        print("No examples found.")
        return 0

    print("Available example configurations:")
    for name, path in examples.items():
        print(f"  {name}: {path}")

    return 0


def validate_config(args):
    """Validate configuration file."""
    try:
        config = load_config(args.config)
        print(f"Configuration file {args.config} is valid JSON.")

        # Basic validation
        if "Problem" not in config:
            print("Warning: Missing 'Problem' section")
        if "Model" not in config:
            print("Warning: Missing 'Model' section")

        print("Basic validation passed.")
        return 0

    except Exception as e:
        print(f"Configuration validation failed: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
