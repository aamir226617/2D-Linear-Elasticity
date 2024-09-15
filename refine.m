function refine(n, degP, dof_Nod)

% Local Nodes on the edges
LocalNodsOnEdge = ElementLocalNumberingScheme(degP);
% number of edges per element (quadrilateral element)
nEdgs=size(LocalNodsOnEdge, 1);

% Connectivity matrix (arranged in local order) 
nodM = [1 2 3 4];                                               % last node is repeated as the first node to represent the edges

% Boundary type
bcTyp = [1 2 3 -1];

% Coordinate matrix (arranged as per node number)
% corM = [0 0;
%            1 0;
%          1/2 2;
%             0 2];
corM = [0 0;   2 0;   2 1;    0 1];

% Constructing connectivity matrix
[nElems, nods_elem]=size( nodM );             % total no of elemets and no of nodes per element 
conM=zeros( nElems , nods_elem * dof_Nod);
for c = 1 : nods_elem                 
    conM( : , [(2*c-1) 2*c]) = [(2*nodM( : , c) - 1)  2*nodM( : , c)];
end
    
save meshInfo1.mat nodM conM corM bcTyp

% number of elements
nElems = size( nodM , 1 );
 
% number of nodes
nNods = size( corM , 1 );
 
% refine p times
if n > 1
 for p = 2 : n
 
 % construting temp_info as empty
 temp_info=[];
 i = 0;                                                                   % counter for temp_info
 
 % we will go to each element 
  for k =1 : nElems
     
     % go to all edges of the element
     for j =1 : nEdgs
         
         y=0;
%          v=nodM(i, j : (j+1));                                        % node numbers on i^th element and j^th edge
         v = nodM( k , LocalNodsOnEdge( j , : ));
         if(~isempty( temp_info ))             
             y = find_modified( v , temp_info );                 % (row no. in temp_info) edge containing nodes in v already refined and the node number is y (if not zero)              
         end
         
         if y == 0                                                         % it does not exist, then calculate new node coords from the edge of parent cell
             nNods = nNods+1;
             nod( j ) = nNods;                                        % used for defining connectivity     
             x = (1/2) * sum( corM( v , : ));                           % coord of the new node as average of end nodes
             corM(nNods,:)=x;
             
             temp_info( i + 1 , : ) = [ nNods  v ];
              i = i + 1;
         else
             nod( j ) = temp_info( y , 1 );                             % if it has already been refined, then first column of temp_info of that row yields node no
         end        
        
     end                                                                    % continue for all edges
     
     % add the mid node
     nNods = nNods + 1;
     nod( j + 1) = nNods;                                           % fifth node in the parent mid point
     x = (1/4) * sum( corM( nodM( k , : ) , : ) );               % average of all  the nodes on the parent cell
     corM( nNods , : ) = x;
     
     % now append connectivity matrix
     nElems = nElems + 1;
     nodM( nElems , : ) = [ nod(1) nodM( k , 2) nod( 2 ) nod( 5 ) ];
     bcTyp( nElems , : ) = [ bcTyp( k , 1 ) bcTyp( k , 2) 0 0];
     nElems = nElems + 1;
     nodM( nElems , : ) = [ nod( 5 ) nod( 2 ) nodM( k , 3 ) nod( 3 )];
     bcTyp( nElems , : ) = [ 0  bcTyp( k , 2 )  bcTyp( k , 3)  0];
     nElems = nElems + 1;
     nodM( nElems , : ) = [nod( 4 ) nod( 5 ) nod( 3 ) nodM( k , 4 ) ];
     bcTyp( nElems , : ) = [0 0 bcTyp( k , 3 ) bcTyp( k , 4 ) ];
     nodM( k , : ) = [nodM( k , 1 ) nod( 1 ) nod( 5 ) nod( 4 ) ];     
     bcTyp( k , : ) = [ bcTyp( k , 1 ) 0 0 bcTyp( k , 4 ) ];
     
%      temp_info
  end
  
% Constructing connectivity matrix
[nElems, nods_elem]=size( nodM )  ;           % total no of elemets and no of nodes per element 
conM=zeros( nElems , nods_elem * dof_Nod) ;
for c = 1 : nods_elem                 
    conM( : , [(2*c-1) 2*c]) = [(2*nodM( : , c) - 1)  2*nodM( : , c)];
end

  fname = sprintf('meshInfo%d.mat', p);
 save(fname, 'nodM', 'conM', 'corM', 'bcTyp')
 end
end
% % TO SHOW MESH
  temp_info=[ ];
 for i = 1 : nElems
     for j = 1 : nEdgs         
         y=0;
%          v=conM(i, j : (j+1));                                        % node numbers on i^th element and j^th edge
         v=nodM(i, LocalNodsOnEdge( j , : ));
         
         if(~isempty(temp_info))             
             y=find_modified(v, temp_info);                 % (row no. in temp_info) edge containing nodes in v already refined and the node number is y (if not zero)              
         end
         if y == 0
             xy=corM(v , :);
             figure(1)
             plot(xy(:,1), xy(: ,2));
             hold on                     
         end         
     end      
     xy=corM(nodM(i , : ), :);
     xc=sum(xy(:,1)) / size(xy, 1);
     yc=sum(xy(:,2)) / size(xy, 1);
     text(xc, yc, num2str(i),'color', 'b');
 end
 for i = 1 : size(corM, 1)
     text(corM(i , 1), corM(i , 2), num2str(i),'color', 'r');
 end
  mytitle = sprintf('RefinementStep-%d', n);
  title(mytitle);
     set(findall(gcf,'-property','FontSize'),'FontSize',20)
  
 
 
 