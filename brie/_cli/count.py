import click

import brie

def print_help(ctx, param, value):
    if value is False:
        return
    click.echo(ctx.get_help())
    ctx.exit()

@click.command()
# @click.pass_context
@click.option("--gffFile", "-a", required=True, default=None,
              type=click.Path(exists=False),
              help="GTF/GFF3 file for gene and transcript annotation")
@click.option("--samList", "-S", required=True, default=None,
              type=click.Path(exists=False),
              help=("A tsv file containing sorted and indexed bam/sam/cram "
                    "files. No header line; file path and cell id (optional)"))
@click.option("--outDir", "-o", required=True, default=None,
              type=click.Path(exists=False),
              help="Full path of output directory [default: $samList/brieCOUNT]")
@click.option("--nproc", "-p", default=1, type=int, show_default=True,
              help="Number of subprocesses")
# @click.option("--add-premRNA", default=False, is_flag=True,
#               help="Add the pre-mRNA as a transcript")
@click.option('--help', '-h', is_flag=True, expose_value=False,
              is_eager=False, callback=print_help, help="Print help message")
def count(GFF, samList, outDir, nproc):
    """Count isoform reads from bam files."""
    add_premRNA = False
    brie.bin.count(GFF, samList, outDir, nproc, add_premRNA)