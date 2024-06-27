import pandas as pd
import matplotlib.pyplot as plt
import sys  # Import the sys module to access command-line arguments

# Check if the user provided a filename as a command-line argument
if len(sys.argv) < 2:
    print("Usage: python3 plot_result.py <filename>")
    sys.exit(1)  # Exit the script if no argument is provided

# The second argument in sys.argv list is the filename (the first one is the script name)
file_name = sys.argv[1]

# Define the full path correctly, removing any misformatted quotes
file_root_path = "/home/ycliao/sphinxsys/build/tests/extra_source_and_tests/test_2d_velocity_poiseuille_flow/bin/"
file_path = file_root_path + file_name

# Load the data from the CSV file
try:
    data = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Error: The file '{file_path}' does not exist.")
    sys.exit(1)

# Ensure the data is in a numpy array format when passing to Matplotlib
position_y = data['Position Y'].to_numpy()  # Convert to numpy array
parabolic_velocity_x = data['Parabolic Velocity X'].to_numpy()
velocity_x = data['Velocity X'].to_numpy()

# Plotting position Y vs parabolic_velocity and velocity
plt.figure(figsize=(12, 6))
plt.plot(position_y, parabolic_velocity_x, 'b-', label='Parabolic Velocity X')
plt.plot(position_y, velocity_x, 'r--', label='Velocity X')
plt.xlabel('Position Y')
plt.ylabel('Velocity X')
plt.title('Position Y vs Velocity X')
plt.grid(True)
plt.legend()
plt.show()
