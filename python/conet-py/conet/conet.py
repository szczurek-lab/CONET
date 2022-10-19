import subprocess
from dataclasses import dataclass
from typing import List, Tuple, Optional

from conet.conet_parameters import CONETParameters
from conet.generative_model.utils import get_logger

logger = get_logger(__name__)


@dataclass
class CONET:
    bin_path: str  # bin_path - path to directory with CoNET binaries
    output_path: Optional[str]
    def infer_tree(self, parameters: CONETParameters):
        if not self.output_path:
            self.output_path = self.bin_path
        try:
            args: List[Tuple[str, str]] = parameters.to_arg_value_pairs()
            cmd = [self.bin_path] + [i for t in args for i in t]
            print(f"Calling CONET executable with args: {args}")
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            with open(f"{self.output_path}CONET_LOG", 'w') as log:
                while process.poll() is None:
                    l = process.stdout.readline()
                    logger.info(f"CONET log: {l}")
                    log.write(f"{l}\n")
        except subprocess.SubprocessError as e:
            logger.error(f"Status : FAIL {e.returncode} {e.output} {e.stdout} {e.stderr}")
