function plot_benchmark(name,datadir,cases,energies,refinement,fixedpoint,cfg,figdir)
%PLOT_BENCHMARK Actual curves, speeds and mesh-refinement diagnostics.
% Four knot counts are four prefixes at fixed H, not four refinement meshes.
% Every SPD coordinate is the actual entry x11,x12,x22. Insets or altered
% position scales are not used in the trajectory comparisons.
if ~exist(figdir,'dir')
    mkdir(figdir);
end
if strcmp(name,'sphere')
    methods = {'reference','interlocking','single_chart','radial_projection'};
else
    methods = {'reference','interlocking','single_chart'};
end
labels={'Reference spline','Fixed-point C^2','Initial-point chart','Radial projection'};
colors=[.1,.1,.1;.05,.35,.75;.83,.2,.16;.04,.55,.28];
styles={'--','-','-.',':'};

% Put all requested knot counts together, with readable panel labels.
% f = new_figure(cfg,[5.8,8]);
% h = [];
for k = 1:numel(cfg.knot_counts)
    N = cfg.knot_counts(k);
    D = load(fullfile(datadir,sprintf('%s_N%d_curves.mat',name,N)));
%     ax = subplot(2,2,k,'Parent',f);
%     ax = axes('Parent',f);
%     hh = draw_curve(ax,name,D,methods,labels,colors,styles);
%     title(ax,sprintf('%d new knots (%d total), T=%g',N,N+1,D.t(end)),'FontSize',9);
%     if k == 1
%         h = hh;
%     end
    % Detailed individual views are supplied in addition to the overview.
    g = new_figure(cfg,[6,5.3]);
    a = axes('Parent',g);
    draw_curve(a,name,D,methods,labels,colors,styles);
%     title(a,sprintf('%s: %d new knots',name,N));
    legend(a,'show','Location','best','FontSize',8);
    save_pdf(g,fullfile(figdir,sprintf('%s_N%d_curve.pdf',name,N)));
    plot_speeds(name,D,methods,labels,colors,styles,cfg,figdir,N);
%     if strcmp(name,'spd2')
%         plot_components(D,methods,labels,colors,styles,cfg,figdir,N);
%     end
end
% legend(h,labels(1:numel(methods)),'Position',[.1,.945,.8,.045], ...
%     'NumColumns',2,'FontSize',8,'Box','off');
% save_pdf(f,fullfile(figdir,[name,'_all_curves.pdf']));
% N = max(cfg.knot_counts);
% copyfile(fullfile(figdir,sprintf('%s_N%d_speeds.pdf',name,N)), ...
%     fullfile(figdir,[name,'_all_speeds.pdf']));
% if strcmp(name,'spd2')
%     copyfile(fullfile(figdir,sprintf('spd2_N%d_components.pdf',N)),fullfile(figdir,'spd2_components.pdf'));
% end

% These errors belong to a FIXED reference, refined within every segment.
f = new_figure(cfg,[6,7.4]);
fields = {'position','velocity','acceleration'};
powers = [5,4,3];
for j = 1:3
    ax = axes('Parent',f,'Position',[.16,.69-.30*(j-1),.79,.22]);
    hold(ax,'on');
    h = [];
    for k = 2:numel(methods)
        selected = refinement(strcmp(refinement.method,methods{k}),:);
        h(k-1) = loglog(ax,selected.N,selected.(fields{j}),[styles{k},'o'], ...
            'Color',colors(k,:),'MarkerSize',4);
    end
    own = refinement(strcmp(refinement.method,'interlocking'),:);
    v = own.(fields{j});
    anchor = max(1,height(own)-1);
    guide = v(anchor)*(own.N/own.N(anchor)).^-powers(j);
    loglog(ax,own.N,guide,'k:');
    set(ax,'XScale','log','YScale','log','XTick',own.N,'FontSize',9);
%     grid(ax,'on');
    xlabel(ax,'Interpolation intervals');
    ylabel(ax,['Maximum ',fields{j},' error']);
    text(ax,.03,.08,sprintf('guide: h^{%d}',powers(j)),'Units','normalized','FontSize',8);
%     if j == 1
%         legend(ax,h,labels(2:numel(methods)),'Position',[.13,.93,.83,.05],'NumColumns',2,'FontSize',8);
%     end
end
save_pdf(f,fullfile(figdir,[name,'_maximum_errors.pdf']));

f = new_figure(cfg,[5.8,3.7]);
ax = axes('Parent',f);
hold(ax,'on');
loglog(ax,fixedpoint.N,fixedpoint.initial_jump,'o--','Color',[.35,.6,.8]);
loglog(ax,fixedpoint.N,fixedpoint.first_jump,'s-','Color',colors(2,:));
loglog(ax,fixedpoint.N,fixedpoint.final_jump,'x:','Color',[.6,.2,.7]);
a = max(1,height(fixedpoint)-1);
nn = fixedpoint.N;
loglog(ax,nn,fixedpoint.initial_jump(a)*(nn/nn(a)).^-3,'k--');
loglog(ax,nn,fixedpoint.first_jump(a)*(nn/nn(a)).^-5,'k:');
set(ax,'XScale','log','YScale','log','XTick',nn);
% grid(ax,'on');
xlabel(ax,'Interpolation intervals');
ylabel(ax,'Maximum actual acceleration jump');
legend(ax,'Initializer','First correction','Returned output','Order 3','Order 5','Location','best','FontSize',8);
save_pdf(f,fullfile(figdir,[name,'_acceleration_jump.pdf']));

f = new_figure(cfg,[5.8,5.9]);
% ax = subplot(2,1,1,'Parent',f);
ax = axes('Parent',f);
hold(ax,'on');
for N = cfg.knot_counts
    D = load(fullfile(datadir,sprintf('%s_N%d_data.mat',name,N)));
    H = D.history;
    semilogy(ax,H(:,1),max(H(:,2),1e-16),'o-','DisplayName',sprintf('%d new knots',N));
end
set(ax,'YScale','log');
% grid(ax,'on');
xlabel(ax,'Fixed-point correction k');
ylabel(ax,'Maximum actual jump');
yline(ax,cfg.fixed_point_tolerance,'k:','tolerance','HandleVisibility','off');
% legend(ax,'show','FontSize',8,'Location','best');
% ax = subplot(2,1,2,'Parent',f);
% valid = isfinite(fixedpoint.local_factor);
% nn = fixedpoint.N(valid);
% qq = fixedpoint.local_factor(valid);
% loglog(ax,nn,qq,'o-','Color',colors(2,:));
% hold(ax,'on');
% loglog(ax,nn,qq(end)*(nn/nn(end)).^-2,'k--');
% grid(ax,'on');
% xlabel(ax,'Interpolation intervals');
% ylabel(ax,'Measured local factor bound');
% legend(ax,'Measured local bound','Second-order guide','Location','best','FontSize',8);
save_pdf(f,fullfile(figdir,[name,'_fixed_point_convergence.pdf']));
% The exact energy values and pairwise separations remain in CSV tables.
assert(height(cases)==numel(cfg.knot_counts) && ~isempty(energies));
end

function h = draw_curve(ax,name,D,methods,labels,colors,styles)
hold(ax,'on');
h = gobjects(1,numel(methods));
if strcmp(name,'sphere')
    [x,y,z] = sphere(16);
    mesh(ax,x,y,z,'FaceColor','none','EdgeColor',[.88,.89,.90], ...
        'EdgeAlpha',.3,'HandleVisibility','off');
end
for k = 1:numel(methods)
    q = components(name,D.([methods{k},'_points']));
    h(k) = plot3(ax,q(:,1),q(:,2),q(:,3),styles{k},'Color',colors(k,:), ...
        'LineWidth',1.5+(k==2)*.3,'DisplayName',labels{k});
end
% N=6 or 10 need not divide the display grid. Mark the ACTUAL IVP knots,
% rather than nearest points on the sampled plotting polyline.
q = components(name,D.knot_points);
scatter3(ax,q(:,1),q(:,2),q(:,3),16,'k','filled','HandleVisibility','off');
scatter3(ax,q(1,1),q(1,2),q(1,3),50,'k','s','HandleVisibility','off');
scatter3(ax,q(end,1),q(end,2),q(end,3),65,'m','*','HandleVisibility','off');
% if strcmp(name,'sphere')
%     axis(ax,'equal');
%     xlabel(ax,'x_1');
%     ylabel(ax,'x_2');
%     zlabel(ax,'x_3');
% else
%     xlabel(ax,'x_{11}');
%     ylabel(ax,'x_{12}');
%     zlabel(ax,'x_{22}');
% end
if strcmp(name,'spd2')
    xlabel(ax,'x_{11}');
    ylabel(ax,'x_{12}');
    zlabel(ax,'x_{22}');
end
view(ax,-56,24);
% grid(ax,'on');
set(ax,'FontSize',8);
end

% function plot_components(D,methods,labels,colors,styles,cfg,figdir,N)
% f = new_figure(cfg,[5.8,6.8]);
% fields = {'x_{11}','x_{12}','x_{22}'};
% for j = 1:3
%     ax = subplot(3,1,j,'Parent',f);
%     hold(ax,'on');
%     for k = 1:numel(methods)
%         q = components('spd2',D.([methods{k},'_points']));
%         plot(ax,D.t,q(:,j),styles{k},'Color',colors(k,:),'LineWidth',1.4,'DisplayName',labels{k});
%     end
%     grid(ax,'on');
%     xlabel(ax,'physical time');
%     ylabel(ax,fields{j});
%     set(ax,'FontSize',9);
%     if j == 1
%         legend(ax,'show','FontSize',8,'Location','best');
%     end
% end
% save_pdf(f,fullfile(figdir,sprintf('spd2_N%d_components.pdf',N)));
% end

function plot_speeds(name,D,methods,labels,colors,styles,cfg,figdir,N)
f = new_figure(cfg,[5.8,5.6]);
% a = subplot(2,1,1,'Parent',f);
a = axes('Parent',f);
% b = subplot(2,1,2,'Parent',f);
g = new_figure(cfg,[5.8,5.6]);
b = axes('Parent',g);
hold(a,'on');
hold(b,'on');
for k = 1:numel(methods)
    speed = D.([methods{k},'_speed']);
    omit = strcmp(name,'sphere') && strcmp(methods{k},'single_chart') && ~cfg.sphere_show_single_chart_speed;
    if ~omit
        plot(a,D.t,speed,styles{k},'Color',colors(k,:),'LineWidth',1.4,'DisplayName',labels{k});
    end
    if k > 1
        semilogy(b,D.t,max(abs(speed-D.reference_speed),1e-16),styles{k}, ...
            'Color',colors(k,:),'LineWidth',1.4,'DisplayName',labels{k});
    end
end
scatter(a,[D.t(1),D.t(end)],D.reference_speed([1,end]),20,'k','filled','HandleVisibility','off');
set(b,'YScale','log');
% grid(a,'on');
% grid(b,'on');
xlabel(a,'Time');
xlabel(b,'Time');
ylabel(a,'Riemannian speed');
ylabel(b,'Absolute speed error');
legend(a,'show','FontSize',8,'Location','best');
legend(b,'show','FontSize',8,'Location','best');
% save_pdf(f,fullfile(figdir,sprintf('%s_N%d_speeds.pdf',name,N)));
save_pdf(f,fullfile(figdir,sprintf('%s_N%d_speeds1.pdf',name,N)));
save_pdf(g,fullfile(figdir,sprintf('%s_N%d_speeds2.pdf',name,N)));
end

function q = components(name,points)
q = zeros(numel(points),3);
for j = 1:numel(points)
    x = points{j};
    if strcmp(name,'sphere')
        q(j,:) = x(:)';
    else
        q(j,:) = [x(1,1),x(1,2),x(2,2)];
    end
end
end

function f = new_figure(cfg,dimensions)
% f = figure('Color','w','Visible',cfg.figure_visible,'Units','inches','Position',[1,1,dimensions]);
f = figure('Color','w','Visible',cfg.figure_visible,'Units','inches');
end

function save_pdf(f,path)
p = get(f,'Position');
set(f,'PaperUnits','inches','PaperSize',p(3:4),'PaperPosition',[0,0,p(3:4)]);
print(f,path,'-dpdf','-painters');
close(f);
end
