

import subprocess


class CoNET:
    # bin_path - path to directory with CoNET binaries
    def __init__(self, bin_path):
        self.bin_path = bin_path

    def infer_tree(self, parameters, input_data):
        input_data.save(parameters.data_dir)
        try:
            cmd = [self.bin_path, parameters.to_string()]
            print(' '.join(cmd))
            cmd_output = subprocess.run(cmd)
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)
        else:
            pass

        # input_data.clear(parameters.data_dir)
        print("TODO")
