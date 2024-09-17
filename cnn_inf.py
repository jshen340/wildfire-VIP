import os
import numpy as np
import tensorflow as tf
from PIL import Image
import random

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
        return float(filename.split('_')[1].split('K')[0])
    except Exception as e:
        print(f"Error extracting temperature from {filename}: {str(e)}")
        return None

# Load the saved model
model = tf.keras.models.load_model('candle_temperature_model.h5')

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
with open('test_set_filenames.txt', 'r') as f:
    test_files = f.read().splitlines()

# Data directory
data_dir = 'CandleFire3000-5000K'

# Function to predict for a specific image
def predict_for_specific_image(image_filename):
    image_path = os.path.join(data_dir, image_filename)
    predicted_temp = predict_temperature(image_path)
    actual_temp = extract_temperature(image_filename)
    
    if predicted_temp is not None and actual_temp is not None:
        print(f"Image: {image_filename}")
        print(f"Predicted temperature: {predicted_temp:.2f}K")
        print(f"Actual temperature: {actual_temp:.2f}K")
        print(f"Absolute error: {abs(predicted_temp - actual_temp):.2f}K")
        print()
    else:
        print(f"Failed to predict temperature for {image_filename}")
        print()

# Function to run inference on multiple selected images
def run_inference_on_selected_images(image_filenames):
    for filename in image_filenames:
        if filename in test_files:
            predict_for_specific_image(filename)
        else:
            print(f"{filename} is not in the test set.")
            print()

# Main execution
if __name__ == "__main__":
    # Inference on a random image in the test set
    # random_test_image = random.choice(test_files)
    # predict_for_specific_image(random_test_image)      
    
    # Inference on several images
    selected_images = [
        "candle3_3087.5K.png",
        "candle46_4341.7K.png",
        "candle381_5858.3K.png",
        "candle361_6441.7K.png",
        "candle338_7112.5K.png",
        "candle185_8395.8K.png",
        "candle231_9737.5K.png",
    ]
    run_inference_on_selected_images(selected_images)
    