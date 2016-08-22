init;
load data_orbits
load matlab
cr3bp = init_CR3BP('EARTH', 'MOON', default);
% We define an unstable manifold                                            
manifold_branch_unstable  = init_manifold_branch(cst.manifold.UNSTABLE,...
                                                 cst.manifold.EXTERIOR);
                                             
%%Plot examples
for k= 6%:size(data_orbits,1)
    % Plot manifold trajectory
    default.plot.manifold_branch= false;

    % Compute minimum distance trajectory
    t0 = 10;
    manifold_branch_unstable = manifold_branch_computation(cr3bp, halo, manifold_branch_unstable,data_orbits(k,1), t0, default, cst, @odezero_3BSOI);
    user.t0=manifold_branch_unstable.te;
    tspan = [user.t0 user.t0+user.tl];
    earth.event.max_events = data_orbits(k,4);
    [~, yearc_bcp, tarc_bcp, yarc_bcp] = ode78_bcp_event(tspan, manifold_branch_unstable.yve(1:6), cr3bp.mu, degtorad(data_orbits(k,2)), cst.sun.ms, cst.sun.as, cst.sun.omegaS, earth.event);
    
    tarc_bcp = [manifold_branch_unstable.ttraj ; tarc_bcp];
    yarc_bcp = [manifold_branch_unstable.ytraj ; yarc_bcp];
    % Plot trajectory
    figure(2*k);
    hold on
    grid on
    hp(2) = plot(yarc_bcp(:,1)*cr3bp.L, yarc_bcp(:,2)*cr3bp.L, 'Color', rgb('dark green'), 'LineWidth', 1.5);
    if(exist('yearc_bcp', 'var'))
        plot(yearc_bcp(:,1)*cr3bp.L, yearc_bcp(:,2)*cr3bp.L, 'o', 'Color', rgb('dark green'), 'MarkerFaceColor', rgb('dark green'), 'MarkerSize', 3);
    end
    %Sphere
    vt = 0:0.01:2*pi;
    tbsoi      = zeros(2, size(vt,2));
    tbsoi(1,:) = cst.env.em3bsoi * cos(vt);
    tbsoi(2,:) = cst.env.em3bsoi * sin(vt);
    %Position
    vl = [1-cr3bp.mu  0  0]*cr3bp.L;
    plot((vl(1) + tbsoi(1,:)),(vl(2) + tbsoi(2,:)), 'k--');
    axis equal


    figure(2*k+1)
    yarc_bcp_B = Geo_to_Helio_changed(yarc_bcp,tarc_bcp,data_orbits(k,2),cr3bp,cst);
    hold on
    grid on
    axis equal
    hp(2) = plot(yarc_bcp_B(:,1)*cr3bp.L*cst.sun.as, yarc_bcp_B(:,2)*cr3bp.L*cst.sun.as, 'Color', rgb('dark green'), 'LineWidth', 1.5);
    plot(yarc_bcp_B(end,1)*cr3bp.L*cst.sun.as, yarc_bcp_B(end,2)*cr3bp.L*cst.sun.as, 'o',  'Color', rgb('dark green'), 'LineWidth', 1.5);
    plot(yarc_bcp_B(1,1)*cr3bp.L*cst.sun.as, yarc_bcp_B(1,2)*cr3bp.L*cst.sun.as, '*',  'Color', rgb('dark green'), 'LineWidth', 1.5);
end