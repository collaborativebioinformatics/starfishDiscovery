import os

def file_paths(input_dir, paths_out, ext):
    """
    doc here
    """
    with open(paths_out, 'w') as fout:
        for f in os.listdir(input_dir):
            full_path = os.path.join(input_dir, f)
            # check for extensions
            if os.path.isfile(full_path) and f.endswith(ext):
                # get the base name without extension
                base_name = os.path.splitext(f)[0]
                fout.write(f"{base_name}\t{full_path}\n")
