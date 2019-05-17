%% Use treestoolbox from
% https://www.treestoolbox.org/manual/index.html
%% Read neurons
if ismac
    cd('/Users/binxu/Holy_Optical_Imaging/Registrate_Atlas/Mitral_cells/')
elseif ispc
    cd('D:\Light-Sheet-Data-Processing-Script\Registrate_Atlas\Mitral_cells\')
end
fn_list = dir();
neuron_list = {}; l=1;
for i = 1:numel(fn_list)
    name = fn_list(i).name;
    if contains(name, '.swc')
        neuron_list{l} = name;
        tree = load_tree(name,'-r');
        trees{l} = tree;
        l=l+1;
    end
end
%% Plot neurons in fish coordinate
figure(3);clf;hold on
for i = 1:l-1
    name = neuron_list{i};
    tree = trees{i};
    G = digraph(tree.dA);
    plot(G,'XData',tree.X,'YData',tree.Y,'ZData',tree.Z,'ShowArrows','off','EdgeAlpha', 0.4, 'LineWidth', 1, ...
                'Marker', 'none')
    soma_id = find(tree.R==2);
    scatter3(tree.X(soma_id), tree.Y(soma_id), tree.Z(soma_id),49,'filled','MarkerFaceAlpha',0.75)
end
axis equal tight
xlabel('X (med-lat R-L)')
ylabel('Y (rost-caud)')
zlabel('Z (vent-dors)')
view([36.6, 11.6])
%%
figure(19);clf;hold on
for i = [1:6]
    name = neuron_list{i};
    tree = trees{i};
    G = digraph(tree.dA);
    plot(G,'XData',tree.X,'YData',tree.Y,'ZData',tree.Z,'ShowArrows','off','EdgeAlpha', 0.6, 'LineWidth', 1.2, ...
                'Marker', 'none')
    soma_id = find(tree.R==2);
    scatter3(tree.X(soma_id), tree.Y(soma_id), tree.Z(soma_id),49,'filled','MarkerFaceAlpha',0.75)
end
axis equal tight
xlabel('X (med-lat R-L)')
ylabel('Y (rost-caud)')
zlabel('Z (vent-dors)')
view([36.6, 11.6])
%% 
G = digraph(tree.dA);
plot(G,'XData',tree.X,'YData',tree.Y,'ZData',tree.Z,'ArrowSize',0)
axis equal
hold on
soma_id = find(tree.R==2);
scatter3(tree.X(soma_id), tree.Y(soma_id), tree.Z(soma_id),49)
