import fire

from dlib.dplotting.plot_projection import plot_projection

def main(
        model_path,
        cif_path,
        map_path
):
    # Produce and save figure for default mu map
    plot_projection(
        model_path,
        cif_path,
        map_path
    )

    # Produce and save figure for pandda2 mu map

    # Produce and save figure for default sigma map

    # Produce and save figure for pandda2 sigma map

    # Produce and save figure for default z map

    # Produce and save figure for pandda2 z map


    ...

if __name__ == "__main__":
    fire.Fire(main)