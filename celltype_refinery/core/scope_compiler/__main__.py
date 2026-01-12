"""
Scope Compiler CLI

Compiles marker map state rule scopes into explicit label lists.

Usage:
    python -m celltype_refinery.core.scope_compiler compile -i INPUT -o OUTPUT
    python -m celltype_refinery.core.scope_compiler validate -i INPUT
    python -m celltype_refinery.core.scope_compiler show -i INPUT
"""

import argparse
import json
import sys
from pathlib import Path

from .compiler import ScopeCompiler
from .config import ScopeCompilerConfig
from .export import format_scope_summary
from .tree import MarkerMapTree
from .validator import validate_marker_map


def cmd_compile(args):
    """Compile scopes and write output."""
    config = ScopeCompilerConfig(output_version=args.version or "v5")
    compiler = ScopeCompiler(config)

    try:
        compiler.compile_file(Path(args.input), Path(args.output))
        print(f"Successfully compiled to {args.output}")
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_validate(args):
    """Validate scope paths without writing."""
    with open(args.input, "r") as f:
        marker_map = json.load(f)

    tree = MarkerMapTree(marker_map)
    errors = validate_marker_map(marker_map, tree)

    if errors:
        print("Validation FAILED:")
        for error in errors:
            print(f"  - {error}")
        return 1
    else:
        print("Validation PASSED: All scope paths are valid")
        return 0


def cmd_show(args):
    """Show resolved scopes (dry run)."""
    with open(args.input, "r") as f:
        marker_map = json.load(f)

    config = ScopeCompilerConfig()
    compiler = ScopeCompiler(config)

    try:
        compiled = compiler.compile(marker_map)
        print(format_scope_summary(compiled))
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Compile marker map state rule scopes into explicit label lists"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # compile command
    compile_parser = subparsers.add_parser(
        "compile", help="Compile scopes and write output JSON"
    )
    compile_parser.add_argument(
        "-i", "--input", required=True, help="Input marker map JSON (v4.1 format)"
    )
    compile_parser.add_argument(
        "-o", "--output", required=True, help="Output marker map JSON (v5 format)"
    )
    compile_parser.add_argument(
        "-v", "--version", default="v5", help="Version string for output (default: v5)"
    )
    compile_parser.set_defaults(func=cmd_compile)

    # validate command
    validate_parser = subparsers.add_parser(
        "validate", help="Validate scope paths without writing"
    )
    validate_parser.add_argument(
        "-i", "--input", required=True, help="Input marker map JSON to validate"
    )
    validate_parser.set_defaults(func=cmd_validate)

    # show command
    show_parser = subparsers.add_parser("show", help="Show resolved scopes (dry run)")
    show_parser.add_argument(
        "-i", "--input", required=True, help="Input marker map JSON"
    )
    show_parser.set_defaults(func=cmd_show)

    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
