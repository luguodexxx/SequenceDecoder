#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/4/10 10:00 PM
__author__ = 'Zhou Ran'

import sys
import click
from .preparedat import SeqCrawler


@click.group()
def cli():
    """Prepare the data to  implement machine/deeping learning. Version: 1.0.0"""
    pass


@click.command()
@click.option('--genome',
              type=str,
              help="The genome fasta file, must be faidxed."
              )
@click.option('--bed',
              type=str,
              help="The ensembl gene id only for now."
              )
@click.option('--fo',
              type=str,
              help="The name of output file."
              )
@click.option('--left',
              type=int,
              help="The left extended length."
              )
@click.option('--right',
              type=int,
              help="The right extended length."
              )
@click.option('--sample',
              type=int,
              default=None,
              help='Sample the input information.')
@click.option('--label',
              type=int,
              default=None,
              help='the label for classification.')
@click.option('--site',
              is_flag=True,
              help='The site mode foe fetch information.')
def bed(genome, bed, fo, left, right, label, sample, site):
    """
    Fetch the sequence from genome
    """
    if not all([genome, bed, fo, left, right]):
        cli(['bed', '--help'])
        sys.exit(1)

    if site:
        file = SeqCrawler.readbed(bed, left, right, label=label, sample=sample, genome=genome)
    else:
        file = SeqCrawler.readbedfromregion(bed, left, right, label=label, sample=sample, genome=genome)

    file.write(fo)


@click.command()
@click.option('--fasta',
              type=str,
              help="The target fasta file."
              )
@click.option('--fo',
              type=str,
              help="The name of output file."
              )
@click.option('--max',
              type=int,
              help="The max length were choosed."
              )
@click.option('--sample',
              type=int,
              help='Sample the input information.')
@click.option('--label',
              type=int,
              help='the label for classification.')
def fasta(fasta, fo, max, sample, label):
    """
    Fetch the sequence from genome
    """

    if not all([fasta, fo, max]):
        cli(['fasta', '--help'])
        sys.exit(1)

    SeqCrawler.randomfasta(fasta, max, sample, label, fileout=fo)


cli.add_command(bed)
cli.add_command(fasta)

if __name__ == '__main__':
    cli()
