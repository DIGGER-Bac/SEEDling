import argparse
import experiment
import logging
logging.basicConfig(level=logging.INFO)



def main(args:argparse.Namespace) -> None:
    logging.info("Starting...")
    experiment.main(args.config)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SEEDling - prediction of seed regions", 
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        "-c", "--config", required=True, type=str, help="Path to the config file (.yml)")

    args = parser.parse_args()

    main(args)