import os
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models
from sklearn.model_selection import train_test_split
from PIL import Image
import random
from tensorflow.keras import regularizers
from tensorflow.keras.optimizers import Adam

# L2 regularization: stabilizes (more consistent results across different trainings), punishes large weights
# Dropout: randomly chooses some inputs to set to 0 during training to prevent overfitting

#CNN for fire shape (temp seems to influence shape more)
#also try fourier transform to only look at color

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
        # Split the filename by underscores and then get the second part
        temp_str = filename.split('_')[1]
        # Remove the ".png" extension and convert to an integer
        return int(temp_str.split('.')[0])
    except Exception as e:
        print(f"Error extracting temperature from {filename}: {str(e)}")
        return None

# Load and preprocess data
data_dir = 'UnitySim5-200'
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
    layers.Conv2D(32, (3, 3), activation='relu', input_shape=(224, 224, 3), 
                  kernel_regularizer=regularizers.l2(0.001)),  # L2 regularization
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu', 
                  kernel_regularizer=regularizers.l2(0.001)),  # L2 regularization
    layers.MaxPooling2D((2, 2)),
    layers.Conv2D(64, (3, 3), activation='relu', kernel_regularizer=regularizers.l2(0.001)),  # L2 regularization
    layers.Flatten(),
    # layers.Dropout(0.5),  # Dropout layer with a dropout rate of 0.5
    layers.Dense(64, activation='relu', kernel_regularizer=regularizers.l2(0.001)),  # L2 regularization
    layers.Dense(1)  # Output layer for regression
])

# Compile the model
model.compile(optimizer=Adam(learning_rate=0.0001), loss='mse', metrics=['mae'])

# Train the model
history = model.fit(X_train, y_train, epochs=100, batch_size=16, validation_split=0.2, verbose=1)
#Using only one picture per temperature, L2, no dropout, LR=0.0001
#35 epochs -- 15.91 MAE
#60 epochs -- 14.12 MAE
#100 epochs -- 12.14 MAE
#150 epochs -- 10.96 MAE !!, 200 epochs -- 10.67 MAE

#Using 2 pictures per temperature, L2, no dropout, LR=0.0001
#35 epochs -- 15.91 MAE
#60 epochs -- 12.14 MAE !!
#100 epochs -- 12.02 MAE
#150 epochs -- 10.53 MAE


# Evaluate the model
test_loss, test_mae = model.evaluate(X_test, y_test, verbose=0)
print(f"Test Mean Absolute Error: {test_mae:.2f}")

# Save the model
model.save('unitySim_temperature_model.h5')

print("Model saved as 'unitySim_temperature_model.h5'")

# Save the test set filenames
with open('unitySim_test_set_filenames.txt', 'w') as f:
    for filename in test_files:
        f.write(f"{filename}\n")

print("Test set filenames saved as 'unitySim_test_set_filenames.txt'")