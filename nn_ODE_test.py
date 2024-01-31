# -*- coding: utf-8 -*-
"""
Deep leaning for differntial equations

Created on Tue Jan 30 11:34:56 2024

@author: user

"""
#========simple 1th order problem================================================
# y = e^{-t/5} sin(t) 
# dy/dt = -1/5*y+cos(t)*e^{-t/5} with y(0)=0 
#
import torch
import numpy as np 
import matplotlib.pyplot as plt 
from torch.autograd import Variable

a = 0
b = 5
N = 100
#---euler method 
y_euler = torch.ones(N)
y_euler[0] = 0
x = torch.Tensor(np.linspace(0,5,100)[:,None])
y=torch.exp(-(x/5))*torch.sin(x) # exact sol 

h = (b-a)/N
for i in range(N-1):
    y_euler[i+1] = y_euler[i] + h* (-1/5.*y_euler[i] +torch.exp(-x[i]/5)*torch.cos(x[i]))
    
fig,ax = plt.subplots()     
ax.plot(x.data.numpy(), y_euler.data.numpy(),'orange',label='euler method')
ax.plot(x.data.numpy(), y.data.numpy(),'g--',label='exact sol')
plt.legend() 

#---NN method 
# assume solution y(x) = A + x*NN(x,theta)
#  this form assigns initial condition y(0)=A 
def sigmoid(x):
    return 1/(1+torch.exp(-x))
def sigmoid_grad(x):
    return sigmoid(x)*(1-sigmoid(x))
def neural_network(x,weights,bias):
    s_z = sigmoid(torch.matmul(x,weights[0])+bias[0])
    return torch.matmul(s_z,weights[1])+bias[1]
def dN_dx(weights,x):
    s_z_grad = sigmoid_grad(torch.matmul(x,weights[0])+bias[0])
    mul = torch.mul(weights[0].T,weights[1])
    return torch.matmul(s_z_grad,mul)

weights = [torch.randn((1,10),requires_grad=True),torch.randn((10,1),requires_grad=True)  ]
bias = [torch.randn(10,requires_grad=True),torch.randn(1,requires_grad=True) ]

A= 0 #initial condition 
Psi_t = lambda x: A+ x *neural_network(x,weights,bias) 
f = lambda x, Psi : torch.exp(-x/5.0)*torch.cos(x) -Psi/5.0 # dydt 

def error(x):    
    # loss = (dy/dt - f(t))**2 
    x.requires_grad = True 
    psi = Psi_t(x) 
    ddN = dN_dx(weights,x) 
    Psi_t_x = neural_network(x,weights,bias) + x*ddN # derivative of assumed y(t)
    return torch.mean((Psi_t_x -f(x,psi))**2) 

epochs = 5000
lr = 0.01 
for i in range(epochs):
    loss = error(x) 
    loss.backward()     
    weights[0].data -= lr*weights[0].grad.data 
    weights[1].data -= lr*weights[1].grad.data 
    bias[0].data -= lr*bias[0].grad.data 
    bias[1].data -= lr*bias[1].grad.data 
    
    weights[0].grad.zero_()
    weights[1].grad.zero_()
    bias[0].grad.zero_()
    bias[1].grad.zero_()
    print("loss: ", loss.item())
    
x = torch.unsqueeze(torch.linspace(0,5,100), dim=1)
y = torch.exp(-(x/5))*torch.sin(x)
psi_trial = Psi_t(x) 
fig,ax = plt.subplots()     
ax.plot(x.data.numpy(), psi_trial.data.numpy(),'orange',label='NN method')
ax.plot(x.data.numpy(), y.data.numpy(),'g--',label='exact sol')
plt.legend() 

#------one can use torch.nn-------------------------------------
print('Note: from here, it is not corect !')
model = torch.nn.Sequential(
          torch.nn.Linear(1,10),
          torch.nn.Sigmoid(),
          torch.nn.Linear(10,1),
          torch.nn.Sigmoid() 
    )

x = Variable(torch.unsqueeze(torch.linspace(0,5,100), dim=1),requires_grad=True)
Psi_nn = lambda x: A+ x *model(x) 
f = lambda x, Psi : torch.exp(-x/5.0)*torch.cos(x) -Psi/5.0 # dydt 
optim = torch.optim.SGD(model.parameters(),lr=0.05)
get_loss = torch.nn.MSELoss(reduction='mean')

for i in range(epochs):
    Psi = Psi_nn(x) 
    #Psi_x = Psi + x*ddN
    # how one can write the derivative over input by using autograd ??? 
    # also, output Psi is not a scalar. thus, one cannot use backward... 
    loss = get_loss(Psi,f(x,Psi))  # this is not correct     
    optim.zero_grad()
    loss.backward()
    optim.step()
    print("loss: ", loss.item())
    
x = torch.unsqueeze(torch.linspace(0,5,100), dim=1)
y = torch.exp(-(x/5))*torch.sin(x)
psi_trial = Psi_nn(x) 
fig,ax = plt.subplots()     
ax.plot(x.data.numpy(), psi_trial.data.numpy(),'orange',label='NN method')
ax.plot(x.data.numpy(), y.data.numpy(),'g--',label='exact sol')
plt.legend()     
    