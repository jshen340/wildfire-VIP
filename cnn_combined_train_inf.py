import os
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models
from sklearn.model_selection import train_test_split
from PIL import Image
import random

# Set random seed for reproducibility
np.random.seed(42)
tf.random.set_seed(42)
random.seed(42)

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

# Load and preprocess data
data_dir = 'CandleFire3000-5000K'
image_files = [f for f in os.listdir(data_dir) if f.endswith('.png')]
temperatures = [extract_temperature(f) for f in image_files]

# Remove any None values (failed to extract temperature)
valid_data = [(f, t) for f, t in zip(image_files, temperatures) if t is not None]
image_files, temperatures = zip(*valid_data)

# Split the data into training and testing sets
train_files, test_files, y_train, y_test = train_test_split(image_files, temperatures, test_size=0.2, random_state=42)

# Load images
X_train = np.array([img for img in (load_and_preprocess_image(os.path.join(data_dir, f)) for f in train_files) if img is not None])
X_test = np.array([img for img in (load_and_preprocess_image(os.path.join(data_dir, f)) for f in test_files) if img is not None])

# Convert y_train and y_test to numpy arrays
y_train = np.array(y_train)
y_test = np.array(y_test)

# Define the CNN model
model = models.Sequential([
    layers.Conv2D(32, (3, 3), activation='relu', input_shape=(224, 224, 3)),
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu'),
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu'),
    layers.Flatten(),
    layers.Dense(64, activation='relu'),
    layers.Dense(1)  # Output layer for regression
])

# Compile the model
model.compile(optimizer='adam', loss='mse', metrics=['mae'])

# Train the model
history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2, verbose=1)

# Evaluate the model
test_loss, test_mae = model.evaluate(X_test, y_test, verbose=0)
print(f"Test Mean Absolute Error: {test_mae:.2f}")

# Function to predict temperature for a new image
def predict_temperature(image_path):
    img = load_and_preprocess_image(image_path)
    if img is not None:
        img = np.expand_dims(img, axis=0)
        predicted_temp = model.predict(img)[0][0]
        return predicted_temp
    else:
        return None

# Example usage with a random test image or manually selected image (should be from test set)
random_test_image = "candle199_8804.2K" #random.choice(test_files)
test_image_path = "candle199_8804.2K.png" #os.path.join(data_dir, random_test_image)
predicted_temperature = predict_temperature(test_image_path)
actual_temperature = extract_temperature(random_test_image)

if predicted_temperature is not None and actual_temperature is not None:
    print(f"Test image: {random_test_image}")
    print(f"Predicted temperature: {predicted_temperature:.2f}K")
    print(f"Actual temperature: {actual_temperature:.2f}K")
    print(f"Absolute error: {abs(predicted_temperature - actual_temperature):.2f}K")
else:
    print("Failed to predict temperature or extract actual temperature.")

# Print out all test set images
print("\nAll test set images:")
for test_file in test_files:
    print(test_file)