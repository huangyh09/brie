import click

import brie

def print_help(ctx, param, value):
    if value is False:
        return
    click.echo(ctx.get_help())
    ctx.exit()

@click.command()
# @click.pass_context
@click.option('--help', '-h', is_flag=True, expose_value=False,
              is_eager=False, callback=print_help, help="Print help message")
def quant():
    """Quantification of isoforms and Detect predictable splicing events with 
    covariant"""
    print("Quant CLI: comming soon.")
    # brie._sh.quant()