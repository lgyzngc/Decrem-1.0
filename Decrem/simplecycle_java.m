function [flag]=simplecycle_java();
flag=0;
% clear java
% javaaddpath('D:\Workspaces\jartest\lib\jgraph-5.13.0.0.jar')
% javaaddpath('D:\Workspaces\jartest\lib\jgrapht-core-0.9.0.jar')
% javaaddpath('D:\Workspaces\jartest\lib\jgrapht-demo-0.9.0.jar')
% javaaddpath('D:\Workspaces\jartest\lib\jgrapht-ext-0.9.0.jar')
% javaaddpath('D:\Workspaces\jartest\lib\jgrapht-ext-0.9.0-uber.jar')
% javaaddpath('D:\Workspaces\jartest\lib\jgraphx-2.0.0.1.jar')
% javaaddpath('D:\Workspaces\jartest\simplecyclesofklength.jar')
h=LinearFBA.FindCycinDirectedGraph;
h.similarMatrix();
flag=1;
% clear java
end
