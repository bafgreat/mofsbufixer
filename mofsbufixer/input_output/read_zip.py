import os
import zipfile
import shutil


def list_files_in_zip(zip_path):
    """
    List all files in the given ZIP archive.

    Parameters:
    zip_path (str): Path to the ZIP file.

    Returns:
    list: A list of file names in the ZIP archive.
    """
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        return zip_ref.namelist()


def search_and_copy_from_zip(mof_names, zip_directory, output_dir):
    """
    Copies the relevant files from all .zip files in the specified directory.

    **Parameters:**
        - mof_names (list): List of MOF names to search for in the .zip files.
        - zip_directory (str): Path to the directory containing .zip files.
        - output_dir (str): Path to the directory where extracted files will be saved.

    **Returns:**
        - str: Path to the directory containing the copied files.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all zip files in the specified directory
    for filename in os.listdir(zip_directory):
        if filename.endswith(".zip"):
            file_path = os.path.join(zip_directory, filename)
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                all_mofs = zip_ref.namelist()
                for mof_name in mof_names:
                    expected_filename = f"Experiment_cif/{mof_name}.cif"
                    if expected_filename in all_mofs:
                        print(f"Copying {expected_filename} from {filename}")

                        # Extract the specific file to a temporary path
                        zip_ref.extract(expected_filename, output_dir)

                        # Construct the path to the extracted file
                        extracted_file_path = os.path.join(
                            output_dir, expected_filename)

                        # Construct the final path for the copied file
                        new_file_path = os.path.join(
                            output_dir, f"{mof_name}.cif")

                        # Copy the extracted file to the new location
                        shutil.copy(extracted_file_path, new_file_path)

                        # Remove the temporary extracted directory if it exists
                        shutil.rmtree(os.path.join(
                            output_dir, "Experiment_cif"), ignore_errors=True)
                        # If found, no need to check other MOF names

    return output_dir
