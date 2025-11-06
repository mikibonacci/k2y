#!/usr/bin/env python3
"""
Command line interface for the k2y package.
Provides tools for working with Koopmans functionals calculations.
"""

import os
import sys
import json
from pathlib import Path

import click

@click.group()
@click.version_option()
def main():
    """k2y: Tools interfacing kcw.x and Yambo codes.
    
    This tool helps prepare and analyze kcw.x data to be used within the Yambo code.
    """
    pass


@main.command('run')
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', type=click.Path(), 
              help='Directory to store calculation results')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
def run_calculation(input_file, output_dir, verbose):
    """Run a Koopmans functional calculation.
    
    INPUT_FILE should be a JSON configuration file for the calculation.
    """
    if verbose:
        click.echo(f"Starting calculation with input file: {input_file}")
    
    try:
        with open(input_file, 'r') as f:
            config = json.load(f)
        
        # Here you would implement the actual calculation logic
        click.echo("Running Koopmans functional calculation...")
        
        # Example placeholder for actual implementation
        from k2y.core import run_koopmans_calculation
        results = run_koopmans_calculation(config, output_dir)
        
        click.echo(f"Calculation completed successfully!")
    except Exception as e:
        click.echo(f"Error during calculation: {str(e)}", err=True)
        sys.exit(1)


@main.command('convert')
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_format', type=click.Choice(['json', 'yaml', 'xml']))
@click.option('--output-file', '-o', type=click.Path(), help='Output file path')
def convert_data(input_file, output_format, output_file):
    """Convert data files between formats.
    
    INPUT_FILE is the file to convert.
    OUTPUT_FORMAT is the desired output format.
    """
    if output_file is None:
        output_file = f"{os.path.splitext(input_file)[0]}.{output_format}"
    
    click.echo(f"Converting {input_file} to {output_format} format...")
    # Implementation of the conversion logic would go here
    click.echo(f"File converted and saved to {output_file}")


@main.command('plot')
@click.argument('data_file', type=click.Path(exists=True))
@click.option('--plot-type', '-p', type=click.Choice(['bands', 'dos', 'orbital']), 
              default='bands', help='Type of plot to generate')
@click.option('--output', '-o', type=click.Path(), help='Save plot to file')
@click.option('--show', is_flag=True, help='Display plot')
def plot_results(data_file, plot_type, output, show):
    """Generate plots from calculation results.
    
    DATA_FILE is the file containing calculation results.
    """
    click.echo(f"Generating {plot_type} plot from {data_file}...")
    
    # Implementation of plotting logic would go here
    # Example placeholder:
    # from k2y.plotting import generate_plot
    # generate_plot(data_file, plot_type, output_file=output, show=show)
    
    if output:
        click.echo(f"Plot saved to {output}")
    

@main.command('info')
@click.argument('file', type=click.Path(exists=True), required=False)
def show_info(file):
    """Display information about a calculation file or the environment.
    
    If FILE is provided, shows details about the file.
    Otherwise, shows information about the k2y environment.
    """
    if file:
        click.echo(f"File information for: {file}")
        # Implementation to show file info would go here
    else:
        click.echo("K2Y Environment Information:")
        click.echo(f"Python version: {sys.version.split()[0]}")
        click.echo(f"K2Y version: {get_version()}")
        # Additional environment information


def get_version():
    """Return the package version."""
    try:
        from importlib.metadata import version
        return version("k2y")
    except:
        return "unknown"


if __name__ == '__main__':
    main()