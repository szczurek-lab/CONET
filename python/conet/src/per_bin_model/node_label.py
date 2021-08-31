
#Label for an event tree in the model
class NodeLabel:
    def __init__(self, start, end, copy_number):
        self.start = start
        self.end = end
        self.copy_number = copy_number
    
    def get_event(self):
        return (self.start, self.end)
    
    def __str__(self):
        return "(" + str(self.start) + "," + str(self.end) + ")"
    
    def overlaps(self, node_label):
        return (self.start >= node_label.start and self.start < node_label.end)\
                or (self.end > node_label.start and self.end < node_label.end)
 
