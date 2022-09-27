import argparse

PHYSICAL_MODEL = None


def main(model_name, solar_model_path):
    print(f"my name is {model_name}")
    print(f"my model is in {solar_model_path}")


def create_model():
    pass


def compute_scattering_image():
    pass


if __name__ == "__main__":
    # Parser initialization
    parser = argparse.ArgumentParser(prog="SPGCoronalModeler",
                                     description="It creates images modelling the Thomson scattering of the solar "
                                                 "corona.",
                                     epilog="Further information at https://aliciavr.github.io/thomsonpy")

    # Arguments creation
    parser.add_argument("-n", "--name", help="Name of the new model")
    parser.add_argument("-m", "--model-path", help="Path to the coronal model to be used.")
    parser.add_argument("echo")

    # Subparsers
    subparsers = parser.add_subparsers()
    # Model creator arguments
    model_creator_parser = subparsers.add_parser(name="Model creator",
                                                 description="Creates a pre-model for the solar image.")
    # Imager arguments
    imager_parser = subparsers.add_parser(name="Imager",
                                          description="Generates solar image.")

    # It reads the arguments from the command line
    args = parser.parse_args()
    print(args.name, args.model_path, args.echo)

    name = "default"
    model = "default"
    echo = "default"

    if args.name:
        name = args.name
    if args.model_path:
        model = args.model_path
    if args.echo:
        echo = args.echo
    main(name, model)
