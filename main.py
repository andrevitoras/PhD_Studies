
if __name__ == '__main__':

    import subprocess
    import os

    requirements_file = "installed_packages.txt"

    # Check if the file exists
    if os.path.exists(requirements_file):
        subprocess.run(["pip", "install", "-r", requirements_file])
        print("Packages installed from requirements file.")
    else:
        print(f"Requirements file '{requirements_file}' not found.")
