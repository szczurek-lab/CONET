from dataclasses import dataclass
import random
from typing import Tuple, List

import numpy as np

Event = Tuple[int, int]
RootEvent: Event = (0, 0)


@dataclass(eq=True, frozen=True)
class NodeLabel:
    """
        Label for nodes of event tree.
        Should be interpreted as:
            cells appended to node with this label have copy_number copy number in segment [start, end)
    """
    start: int
    end: int
    copy_number: int

    def get_event(self) -> Event:
        return self.start, self.end

    def __str__(self):
        return "(" + str(self.start) + "," + str(self.end) + ")"

    def overlaps(self, node_label: 'NodeLabel') -> bool:
        return (node_label.start <= self.start < node_label.end) or (node_label.start < self.end < node_label.end)


class EventSampler:
    max_event_length: int = 100

    @staticmethod
    def sample_events(no_loci: int, number_of_events: int) -> List[Event]:
        """
            Returns sample of @number_of_events distinct events for chromosome
            with @no_loci loci
        """
        candidates = [(a, b) for a in range(0, no_loci) for b in range(0, no_loci) if a < b and b - a < EventSampler.max_event_length]

        result = []
        for n in range(0, number_of_events):
            event = EventSampler.__sample_event(candidates)
            candidates.remove(event)
            result.append(event)
        return result

    @staticmethod
    def __sample_event(candidate_events: List[Event]):
        norming_constant = sum([np.exp(-0.05 * (z[1] - z[0])) for z in candidate_events])
        while True:
            label = random.sample(candidate_events, 1)[0]
            if np.random.uniform() < np.exp(-0.05 * (label[1] - label[0])) / norming_constant:
                return label
