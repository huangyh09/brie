
from ..version import __version__

def main():
    print("Welcome to BRIE v%s! Command lines available:\n " %(__version__))
    print("brie-count\n    Counting reads for exon-skipping events from "
          "per cell bam files")
    print("brie-quant\n    Quantify splicing ratio and detecting differential "
          "splicing\n")


if __name__ == "__main__":
    main()
