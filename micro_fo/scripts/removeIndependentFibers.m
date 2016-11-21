function [nodes, fibers] = removeIndependentFibers(nodes,fibers)
%
% Description:
%   Removes any fibers that aren't attached to at least one other
%
%   Changed to be a function that takes and returns fibers and nodes instead
%   of reading and writing from file
%
% By: Dan Fovargue
% Last update: 4/2015
%

if size(nodes,2) ~= 4
    error('Nodes structure wrong size in removeIndependentFibers')
end
if size(fibers,2) ~= 3
    error('Fibers structure wrong size in removeIndependentFibers')
end

% Get each nodes connection count
connect = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
    sum1 = sum(fibers(:,2)==nodes(i,1));
    sum2 = sum(fibers(:,3)==nodes(i,1));
    connect(i) = sum1 + sum2;
end

fiberRemove = zeros(size(fibers,1),1);
nodeRemove = zeros(size(nodes,1),1);

% If both nodes have 1 connection then remove
for i=1:size(nodes,1)
    if(connect(i)==1)
        fiber = find((fibers(:,2)==i));
        if(fiber)
            otherNode = fibers(fiber,3);
            if(connect(otherNode)==1)
                fiberRemove(fiber) = 1;
                nodeRemove(i) = 1;
                nodeRemove(otherNode) = 1;
            end
        end
        fiber = find((fibers(:,3)==i));
        if(fiber)
            otherNode = fibers(fiber,2);
            if(connect(otherNode)==1)
                fiberRemove(fiber) = 1;
                nodeRemove(i) = 1;
                nodeRemove(otherNode) = 1;
            end
        end
    end
end

nodes = nodes(~nodeRemove,:);
fibers = fibers(~fiberRemove,:);

% have to change node indices
for i=1:size(nodes,1)
    
    if(nodes(i,1) ~= i)
        fibers = [fibers(:,1) ...
                  (fibers(:,2)==nodes(i,1))*i + (fibers(:,2)~=nodes(i,1)).*fibers(:,2) ...
                  (fibers(:,3)==nodes(i,1))*i + (fibers(:,3)~=nodes(i,1)).*fibers(:,3)];
        nodes(i,1) = i;
    end
end

end

