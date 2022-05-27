import subprocess
from dataclasses import dataclass

from conet.conet_parameters import CONETParameters
from conet.generative_model.utils import get_logger

logger = get_logger(__name__)


@dataclass
class CONET:
    bin_path: str  # bin_path - path to directory with CoNET binaries

    def infer_tree(self, parameters: CONETParameters):
        try:
            cmd = [self.bin_path] + parameters.to_string()
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            with open(f"{self.bin_path}_LOG", 'w') as log:
                while process.poll() is None:
                    l = process.stdout.readline()
                    logger.info(f"CONET log: {l}")
                    log.write(f"{l}\n")
        except subprocess.SubprocessError as e:
            logger.error(f"Status : FAIL {e.returncode} {e.output} {e.stdout} {e.stderr}")
