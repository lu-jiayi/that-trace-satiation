import pandas as pd
import json
import os

# Ask for the stimuli Excel file name
file_name = input("Enter the stimuli file name (e.g., exp3_stim.xlsx): ")

# Check if the file exists
if not os.path.isfile(file_name):
    print(f"Error: File '{file_name}' not found in the current folder.")
else:
    # Read the file
    excel_data = pd.read_excel(file_name, sheet_name=0)

    # Convert to JSON
    json_data = excel_data.to_json(orient='records', indent=4)

    # Write the JSON data to a new file
    with open("stimuli.json", "w", encoding='utf-8') as json_file:
        json.dump(json.loads(json_data), json_file, ensure_ascii=False, indent=4)

    # Print confirmation
    print("JSON file created successfully: stimuli.json")