

import subprocess


class CoNET:
    # bin_path - path to directory with CoNET binaries
    def __init__(self, bin_path):
        self.bin_path = bin_path

    def infer_tree(self, parameters, input_data):
        input_data.save(parameters.data_dir)
        try:
            cmd = [self.bin_path] + parameters.to_string()
            print(' '.join(cmd))
            process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            while True: 
                output = process.stdout.readline()
                if not output: 
                    break 
                else: 
                    print(output.rstrip()) 
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            pass
            
    def infer_tree(self, parameters):
        try:
            cmd = [self.bin_path] + parameters.to_string()
            print(' '.join(cmd))
            process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            while True: 
                output = process.stdout.readline()
                if not output: 
                    break 
                else: 
                    print(output.rstrip()) 
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            pass

