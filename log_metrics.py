""" Function to log operating metrics (mostly from the collector) on REBS runs. """
import os
import pandas as pd
import re


def log_collector(filename="Metrics.xml"):
    collector_log_dict = {}
    countdown = 0
    shutdown = False
    try:
        with open(filename) as file:
            foldername = (re.findall(r'[/\\][0-9][0-9][0-9][0-9].+[/\\]', os.path.abspath(filename)))[-1].strip('/\\')
            print(f"Folder name: {foldername}")
            for line in file:
                if 'Info Type="Collector"' in line:
                    countdown = 12
                    shutdown = True
                    continue
                if countdown:
                    MetricName = re.search(r'["].*["]', line).group().strip('"')
                    MetricData = re.search(r'[>].*[<]', line).group().strip('><')
                    if MetricData:
                        collector_log_dict[MetricName] = MetricData
                    else:
                        print(f"The metric data on line {13-countdown} after the beginning of the Collector block could not be found. Cancelling operation. Please retry.")
                        return None
                    countdown -= 1
                elif (not countdown) and shutdown:
                    collector_log_dict['Data File'] = foldername
                    return collector_log_dict
    except FileNotFoundError:
        print(f"No file with that name (default is {filename}) was found. Please try again.")


def save_collector(data, master_csv_path='D:/Documents/Huffman Lab/Instruments/REBS/Data/masterCollectorLogs.csv'):
    df = pd.DataFrame(data, index=[0])
    if os.path.exists(master_csv_path):
        df.to_csv(master_csv_path, mode='a', index=False, header=False)
    else:
        df.to_csv(master_csv_path, index=False)
    return "Collector metrics saved!"


def collect_all_data(data_dir='D:\\Documents\\Huffman Lab\\Instruments\\REBS\\Data', filename='Metrics.xml'):
    original_cwd = os.getcwd()
    for base_path, dir_names, file_names in os.walk(data_dir):
        if filename in file_names:
            os.chdir(base_path)
            collector_data = log_collector(filename)
            save_collector(collector_data)
    os.chdir(original_cwd)


if __name__ == "__main__":
    collector_data = log_collector()
    save_collector(collector_data)
