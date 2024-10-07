"""Console script for microorganism_snps."""
import microorganism_snps

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for microorganism_snps."""
    console.print("Replace this message by putting your code into "
               "microorganism_snps.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    


if __name__ == "__main__":
    app()
