

import subprocess


class CONET:
    # bin_path - path to directory with CoNET binaries
    def __init__(self, bin_path):
        self.bin_path = bin_path

    def infer_tree(self, parameters, input_data=None):
        if input_data is not None:
            input_data.save(parameters.data_dir)
        try:
            cmd = [self.bin_path] + parameters.to_string()
            print(' '.join(cmd))
            process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            while process.poll() is None:
                l = process.stdout.readline()
                print(l)
            print(process.stdout.read())
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            pass
            
