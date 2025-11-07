import numpy as np
import os
import gzip

# Create data directory if it doesn't exist
os.makedirs('data', exist_ok=True)

def read_stars_data_fixed(filename):
    """
    Read stars.dat file handling split lines where the 6th column 
    appears on the next line
    """
    # Check if file is gzipped
    if filename.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    data_list = []
    incomplete_row = []
    
    with opener(filename, mode) as f:
        for i, line in enumerate(f):
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split the line and convert to float
            values = line.strip().split()
            
            # If we have an incomplete row from previous line
            if incomplete_row:
                # Add the current values to complete the row
                incomplete_row.extend(values)
                if len(incomplete_row) >= 6:
                    try:
                        row = [float(v) for v in incomplete_row[:6]]
                        data_list.append(row)
                        # If there are extra values, save them for next row
                        if len(incomplete_row) > 6:
                            incomplete_row = incomplete_row[6:]
                        else:
                            incomplete_row = []
                    except ValueError as e:
                        print(f"Error parsing combined line at {i+1}: {e}")
                        incomplete_row = []
                continue
            
            # Check if we have exactly 6 columns
            if len(values) == 6:
                try:
                    row = [float(v) for v in values]
                    data_list.append(row)
                except ValueError as e:
                    print(f"Error parsing line {i+1}: {e}")
            
            # If we have 5 columns, save them and expect the 6th on next line
            elif len(values) == 5:
                incomplete_row = values
            
            # If we have 1 column and no incomplete row, this might be an error
            elif len(values) == 1 and not incomplete_row:
                print(f"Warning: Orphaned single value at line {i+1}: {values[0]}")
            
            # For any other number of columns
            else:
                if len(values) < 6:
                    incomplete_row = values
                else:
                    # More than 6 values, take first 6
                    try:
                        row = [float(v) for v in values[:6]]
                        data_list.append(row)
                        if len(values) > 6:
                            incomplete_row = values[6:]
                    except ValueError as e:
                        print(f"Error parsing line {i+1}: {e}")
    
    return np.array(data_list)

# Read the file
try:
    if os.path.exists('stars.dat'):
        data = read_stars_data_fixed('stars.dat')
    elif os.path.exists('stars.dat.gz'):
        data = read_stars_data_fixed('stars.dat.gz')
    else:
        raise FileNotFoundError("Neither 'stars.dat' nor 'stars.dat.gz' found")
    
    print(f"Successfully loaded {len(data)} rows of data")
    print(f"Data shape: {data.shape}")
    
    # Verify the data looks correct
    print("\nFirst 5 rows:")
    print(data[:5])
    print("\nLast 5 rows:")
    print(data[-5:])
    
except Exception as e:
    print(f"Error reading file: {e}")
    raise

# Save the cleaned data for future use
np.savetxt('data/stars_cleaned.dat', data, fmt='%15.6E')
print(f"\nCleaned data saved to 'data/stars_cleaned.dat'")
