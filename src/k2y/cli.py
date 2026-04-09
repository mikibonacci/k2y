#!/usr/bin/env python3
"""
Command line interface for the k2y package.
Provides tools for converting Koopmans eigenvalues into Yambo QP databases.
"""

import sys
from pathlib import Path

import click


@click.group()
@click.version_option()
def main():
    """k2y: Convert Koopmans (kcw.x) eigenvalues into Yambo QP databases.

    Supports two modes:

    \b
    File mode  – provide ns.db1, a .kho eigenvalue file, and a pw.x input used for the nscf kcw run.
    AiiDA mode – provide Yambo and KCW node PKs from an AiiDA database.
    """
    pass


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------

@main.command('generate')
# --- file-mode inputs -------------------------------------------------------
@click.option('--ns-db1', type=click.Path(), default=None,
              help='Path to Yambo ns.db1 file (SAVE/ directory). '
                   'Required in file mode.')
@click.option('--eval', 'koopmans_eval', type=click.Path(), default=None,
              help='Path to the kcw.x output file containing Koopmans eigenvalues '
                   '(must be readable by ase_koopmans.io, typically a .kho file). '
                   'Required in file mode (or hybrid AiiDA/file mode).')
@click.option('--pwinput', type=click.Path(), default=None,
              help='Path to the pw.x input file used for the KCW calculation. '
                   'Must contain the same k-point grid as the KCW run. Required in file mode.')
# --- AiiDA-mode inputs ------------------------------------------------------
@click.option('--yambo-pk', type=int, default=None,
              help='AiiDA node PK of the Yambo calculation. '
                   'Activates AiiDA mode.')
@click.option('--kcw-pk', type=int, default=None,
              help='AiiDA node PK of the KCW calculation (optional in AiiDA mode). '
                   'If omitted, --eval and --pwinput must be provided.')
# --- common options ---------------------------------------------------------
@click.option('--template-qp', type=click.Path(), default=None,
              help='Path to a custom template ndb.QP file. '
                   'If omitted, the bundled template is used.')
@click.option('--spin', is_flag=True, default=False,
              help='Use the spin-polarised bundled template '
                   '(ignored when --template-qp is given).')
@click.option('--output', '-o', default='ndb.QP',
              help='Output QP database filename. [default: ndb.QP]')
@click.option('--time-rev/--no-time-rev', default=True,
              help='Use time-reversal symmetry when matching k-points. '
                   '[default: --time-rev]')
@click.option('--brute-force/--no-brute-force', default=True,
              help='Fall back to brute-force |k| matching. '
                   '[default: --brute-force]')
# --- verification options ---------------------------------------------------
@click.option('--verify-k', type=int, default=None,
              help='1-based k-point index used for mapping verification.')
@click.option('--verify-tv', type=int, default=None,
              help='1-based top-valence band index used for mapping verification.')
# --- AiiDA store option -----------------------------------------------------
@click.option('--store', is_flag=True, default=False,
              help='Store the output as an AiiDA SinglefileData node '
                   '(AiiDA mode only).')
def generate(ns_db1, koopmans_eval, pwinput,
             yambo_pk, kcw_pk,
             template_qp, spin, output,
             time_rev, brute_force,
             verify_k, verify_tv,
             store):
    """Generate a Yambo QP database from Koopmans eigenvalues.

    \b
    File mode example:
      k2y generate --ns-db1 SAVE/ns.db1 --eval kc.kho --pwinput nscf.in

    \b
    AiiDA mode example (full):
      k2y generate --yambo-pk 1234 --kcw-pk 5678

    \b
    Hybrid mode (AiiDA Yambo + local KCW files):
      k2y generate --yambo-pk 1234 --eval kc.kho --pwinput nscf.in
    """
    from k2y.k2y import KcwQpDatabaseGenerator

    aiida_mode = yambo_pk is not None

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------
    if not aiida_mode:
        missing = []
        if not ns_db1:
            missing.append('--ns-db1')
        if not koopmans_eval:
            missing.append('--eval')
        if not pwinput:
            missing.append('--pwinput')
        if missing:
            raise click.UsageError(
                f"File mode requires: {', '.join(missing)}. "
                "Alternatively, use --yambo-pk to switch to AiiDA mode."
            )
    else:
        # AiiDA mode: if kcw-pk is absent we need the local files
        if kcw_pk is None:
            missing = []
            if not koopmans_eval:
                missing.append('--eval')
            if not pwinput:
                missing.append('--pwinput')
            if missing:
                raise click.UsageError(
                    f"When --kcw-pk is not provided, {', '.join(missing)} are required "
                    "(hybrid mode: AiiDA Yambo + local KCW files)."
                )

    # ------------------------------------------------------------------
    # Build converter
    # ------------------------------------------------------------------
    try:
        if aiida_mode:
            click.echo(f"AiiDA mode: yambo_pk={yambo_pk}, kcw_pk={kcw_pk}")
            converter = KcwQpDatabaseGenerator.from_aiida(
                yambo_node_pk=yambo_pk,
                kcw_node_pk=kcw_pk,
                template_QP_path=template_qp,
                spin=spin,
            )
        else:
            click.echo("File mode")
            converter = KcwQpDatabaseGenerator(
                ns_db1=ns_db1,
                template_QP_path=template_qp,
                spin=spin,
            )

        # Set eigenvalues from file when not provided by AiiDA
        if koopmans_eval is not None:
            click.echo(f"Loading Koopmans eigenvalues from: {koopmans_eval}")
            converter.set_koopmans_eval(path=koopmans_eval)

        # Set k-points from pw.x input when not provided by AiiDA
        if pwinput is not None:
            click.echo(f"Loading k-points from: {pwinput}")
            converter.set_kpoints_from_pwinput(pwinput)

        click.echo(converter.summary())

        # ------------------------------------------------------------------
        # Core workflow
        # ------------------------------------------------------------------
        click.echo("Generating mappings...")
        converter.generate_mappings(time_rev=time_rev, brute_force=brute_force)

        if verify_k is not None and verify_tv is not None:
            click.echo(f"Verifying mappings at k={verify_k}, top_valence={verify_tv}...")
            converter.verify_mappings(k_index=verify_k, top_valence=verify_tv)

        click.echo(f"Writing QP database to: {output}")
        converter.generate_QP_db(output_filename=output)

        if store:
            if not aiida_mode:
                click.echo("Warning: --store is only meaningful in AiiDA mode.", err=True)
            else:
                click.echo("Storing output as AiiDA SinglefileData...")
                node = KcwQpDatabaseGenerator.generate_SinglefileData_from_file(output)
                node.store()
                click.echo(f"Stored with pk={node.pk}")

        click.echo("Done.")

    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# kpoints
# ---------------------------------------------------------------------------

@main.command('kpoints')
@click.option('--ns-db1', type=click.Path(exists=True), required=True,
              help='Path to Yambo ns.db1 file.')
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Write K_POINTS card to this file (default: print to stdout).')
@click.option('--coordinates', type=click.Choice(['crystal', 'tpiba']),
              default='crystal', show_default=True,
              help='Coordinate system for k-points.')
@click.option('--full-bz', is_flag=True, default=False,
              help='Expand irreducible BZ to full Brillouin zone.')
def kpoints(ns_db1, output, coordinates, full_bz):
    """Print the K_POINTS card from a Yambo ns.db1 file.

    Useful for setting up a kcw.x interpolation run that covers
    exactly the same k-point grid as the Yambo calculation.

    \b
    Example:
      k2y kpoints --ns-db1 SAVE/ns.db1 --output kpoints.txt
    """
    from k2y.k2y import KcwQpDatabaseGenerator

    try:
        converter = KcwQpDatabaseGenerator(ns_db1=ns_db1)
        converter.produce_kpoints_for_interpolation(
            filename=output,
            coordinates=coordinates,
            full_BZ=full_bz,
        )
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@main.command('info')
def show_info():
    """Display k2y environment information."""
    try:
        from importlib.metadata import version as pkg_version
        v = pkg_version("k2y")
    except Exception:
        v = "unknown"

    click.echo(f"k2y version : {v}")
    click.echo(f"Python      : {sys.version.split()[0]}")


if __name__ == '__main__':
    main()
