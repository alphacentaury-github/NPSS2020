# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 23:51:38 2024

@author: young
"""

import torch 
import torch.nn as nn 
import torch.optim as optim 
from torchvision import datasets, transforms 
from sklearn.metrics import classification_report 

class MyModule(nn.Module): 
	def __init__(self, num_inputs, num_outputs, hidden_size): 
		super(MyModule, self).__init__() 
		self.linear1 = nn.Linear(num_inputs, hidden_size) 
		self.linear2 = nn.Linear(hidden_size, num_outputs) 

	def forward(self, input): 
		lin = self.linear1(input) 
		output = nn.functional.relu(lin) 
		pred = self.linear2(output) 
		return pred 

# Instantiate the custom module 
my_module = MyModule(num_inputs=28*28, num_outputs=10, hidden_size=20) 

# Define the loss function and optimizer 
criterion = nn.CrossEntropyLoss() 
optimizer = optim.SGD(my_module.parameters(), lr=0.01) 

# Define the transformations for the dataset 
transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5,), (0.5,))]) 

# Load the MNIST dataset 
train_dataset = datasets.MNIST(root='./data', train=True, download=True, transform=transform) 
test_dataset = datasets.MNIST(root='./data', train=False, download=True, transform=transform) 

# Define the data loader 
train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=64, shuffle=True) 
test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=64, shuffle=False) 

# Train the model 
for epoch in range(10): 
	for i, (images, labels) in enumerate(train_loader): 
		images = images.view(-1, 28*28) 
		optimizer.zero_grad() 
		output = my_module(images) 
		loss = criterion(output, labels) 
		loss.backward() 
		optimizer.step() 
	print('Epoch -->',epoch,'-->',loss) 

	

#Test the model 
with torch.no_grad(): 
	y_true = [] 
	y_pred = [] 
	correct = 0
	total = 0
	for images, labels in test_loader: 
		images = images.view(-1, 28*28) 
		output = my_module(images) 
		_, predicted = torch.max(output.data, 1) 
		total += labels.size(0) 
		correct += (predicted == labels).sum() 
		y_true += labels.tolist() 
		y_pred += predicted.tolist() 

	# Accuracy 
	print('Accuracy: {} %'.format(100 * correct / total)) 
	
	# Classification Report 
	report = classification_report(y_true, y_pred) 
	print(report)
