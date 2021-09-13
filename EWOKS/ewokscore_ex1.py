from ewokscore import Task
from ewokscore import load_graph

# Implement a workflow task
class SumTask(
    Task, input_names=["a"], optional_input_names=["b"], output_names=["result"]
):
    def run(self):
        result = self.inputs.a
        if self.inputs.b:
            result += self.inputs.b
        self.outputs.result = result

# Define a workflow
nodes = [
    {"id": "task1", "class": "SumTask", "inputs": {"a": 1}},
    {"id": "task2", "class": "SumTask", "inputs": {"b": 2}},
    {"id": "task3", "class": "SumTask", "inputs": {"b": 1}},
]
links = [
    {"source": "task1", "target": "task2", "arguments": {"a": "result"}},
    {"source": "task2", "target": "task3", "arguments": {"a": "result"}},
]
workflow = {"nodes": nodes, "links": links}

# Execute a workflow (use a proper Ewoks task scheduler in production)
graph = load_graph(workflow)
print(graph)
# varinfo = {"root_uri": "/tmp/myresults"}  # optional
tasks = graph.execute() #varinfo=varinfo)
# print(tasks)
# print("task 1, output_values = ", tasks["task1"].output_values)
# print("task 2, output_values = ", tasks["task2"].output_values)
# print("task 3, output_values = ", tasks["task3"].output_values)