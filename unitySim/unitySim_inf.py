import os
import numpy as np
import tensorflow as tf
from PIL import Image
import random
import matplotlib.pyplot as plt  # Import matplotlib for plotting
import re

# Function to load and preprocess images
def load_and_preprocess_image(file_path):
    try:
        with Image.open(file_path) as img:
            img = img.convert('RGB')  # Convert to RGB mode
            img = img.resize((224, 224))
            img_array = np.array(img)
            img_array = img_array / 255.0  # Normalize pixel values
        return img_array
    except Exception as e:
        print(f"Error processing image {file_path}: {str(e)}")
        return None

# Function to extract temperature from filename
def extract_temperature(filename):
    try:
        # Split the filename by underscores and then get the second part
        temp_str = filename.split('_')[1]
        # Remove the ".png" extension and convert to an integer
        return int(temp_str.split('.')[0])
    except Exception as e:
        print(f"Error extracting temperature from {filename}: {str(e)}")
        return None

# Load the saved model
model = tf.keras.models.load_model('unitySim_temperature_model.h5') 

# Function to predict temperature for a new image
def predict_temperature(image_path):
    img = load_and_preprocess_image(image_path)
    if img is not None:
        img = np.expand_dims(img, axis=0)
        predicted_temp = model.predict(img)[0][0]
        return predicted_temp
    else:
        return None

# Load test set filenames
with open('unitySim_test_set_filenames.txt', 'r') as f:
    test_files = f.read().splitlines()

# Data directory
data_dir = 'UnitySim5-200'

# Function to predict for a specific image
def predict_for_specific_image(image_filename):
    image_path = os.path.join(data_dir, image_filename)
    predicted_temp = predict_temperature(image_path)
    actual_temp = extract_temperature(image_filename)
    
    if predicted_temp is not None and actual_temp is not None:
        absolute_error = abs(predicted_temp - actual_temp)
        return image_filename, predicted_temp, actual_temp, absolute_error
    else:
        return None, None, None, None

# Function to run inference on multiple selected images
def run_inference_on_selected_images(image_filenames):
    errors = []
    actual_temps = []
    predicted_temps = []
    valid_filenames = []  # To keep track of valid filenames
    
    for filename in image_filenames:
        if filename in test_files:
            result = predict_for_specific_image(filename)
            if result[0] is not None:
                image_filename, predicted_temp, actual_temp, absolute_error = result
                errors.append(absolute_error)
                actual_temps.append(actual_temp)
                predicted_temps.append(predicted_temp)
                valid_filenames.append(image_filename)  # Store the valid filename
                print(f"Image: {image_filename}")
                print(f"Predicted temperature: {predicted_temp:.2f}")
                print(f"Actual temperature: {actual_temp:.2f}")
                print(f"Absolute error: {absolute_error:.2f}")
                print()
            else:
                print(f"Failed to predict temperature for {filename}")
                print()
        else:
            print(f"{filename} is not in the test set.")
            print()
    
    return valid_filenames, errors, actual_temps, predicted_temps


def extract_number(filename):
    """Extracts the numeric part of the filename, ignoring any digits that follow 'Temp' or 'Temp1'."""
    # Use regex to match the pattern after 'Temp' or 'Temp1' and extract only the temperature digits
    match = re.search(r'Temp1?_(\d+)', filename)  # Match 'Temp_' or 'Temp1_' followed by digits
    return int(match.group(1)) if match else float('inf')  # Return the extracted number or infinity if not found

def sort_filenames_from_file(file_path):
    """Reads filenames from a file and returns them sorted by the numeric part."""
    try:
        with open(file_path, 'r') as f:
            file_names = f.read().splitlines()  # Read lines and strip newline characters
        
        # Sort the filenames based on the extracted number
        sorted_file_names = sorted(file_names, key=extract_number)

        return sorted_file_names
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return []

# Main execution
if __name__ == "__main__":
    # Inference on several images
    # selected_images = [
    #     "Temp_006.png",
    #     "Temp_007.png",
        
    #     "Temp1_012.png",
    #     "Temp1_017.png",
        
    #     "Temp1_022.png",
    #     "Temp_026.png",
    # ]
    
    #inference on all of the test set, sorted by temperature, then plotted
    
    sorted_filenames = sort_filenames_from_file('unitySim_test_set_filenames.txt')
    print(sorted_filenames)
    
    valid_filenames, errors, actual_temps, predicted_temps = run_inference_on_selected_images(sorted_filenames)

    # Plotting the absolute error
    plt.figure(figsize=(10, 6))
    plt.bar(valid_filenames, errors, color='orange')
    plt.xlabel('Image Filenames (ascending temperature)')
    plt.ylabel('Absolute Error')
    plt.title('Absolute Error of Temperature Predictions')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
