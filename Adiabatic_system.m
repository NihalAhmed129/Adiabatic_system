clear
L = [50 50]; % dimensions of the box
N = 1000; % number of particles
D = 1/4; % diameter of particles
r = rand(N,2)*(L(1,1)-D/2); % define initial positions
v = randn(N,2); % define initial velocities
dt = 0.1; % time step

T = 10000; % number of steps
Q = 100; % ensemble over time

M = 10;
d_vx = 1; min_vx = -5;
d_vy = 1; min_vy = -5;
d_rx = 2*L(1,1)/(M - 1); min_rx = 0;
d_ry = 2*L(1,2)/(M - 1); min_ry = 0;
d_ps = d_vx*d_vy*d_rx*d_ry;
n_ps(Q,M^4) = 0;
entropy(T,1) = 0; % entropy calculated from phase space

vL=[0.01 0];

for t=1:T
    
    d_rx = L(1,1)/9;
    d_ry = L(1,2)/9;
    d_vrxy = d_vx*d_vy*d_rx*d_ry;
    
    pressure_q(1:Q,1) = 0;
    temperature_q(1:Q,1) = 0;
    
    for q=1:Q
          % This part defines the collisions between balls
    [sx,nx] = sort(r(:,1)); % sort the x coordinate
    for kn = 1 : N-1
        n = nx(kn); % scan particles from left to right
        for km = kn+1:N
            m = nx(km);
            if (sx(km) - sx(kn) > D) % stop if x's differ more than D
                break
            else
                rp = r(n,:)-r(m,:); % vector between two balls
                nrp = norm(rp);
                if (nrp < D) % make sure that particles are close enough
                    rv = v(n,:)-v(m,:); % velocity difference vector
                    if (rv*rp'<0) % make sure that particles are moving towards each other
                        v(n,:) = v(n,:) - (rv*rp')*rp/nrp^2; % bounce
                        v(m,:) = v(m,:) + (rv*rp')*rp/nrp^2; % bounce
                    end
                end
            end
        end
    end
    
    %this part calculates the collusions of particles with wall and
    %pressure
     nv = (v(:,1) > 0).*(r(:,1) > L(1,1)-D/2); v(nv==1,1) = -v(nv==1,1)+2*vL(1,1); pressure_q(q,1)=pressure_q(q,1)-2*sum(v(nv==1,1))/dt/L(1,2)-2*sum(nv==1)*vL(1,1)/dt/L(1,2);
     nv = (v(:,2) > 0).*(r(:,2) > L(1,2)-D/2); v(nv==1,2) = -v(nv==1,2)+2*vL(1,2); pressure_q(q,1)=pressure_q(q,1)-2*sum(v(nv==1,2))/dt/L(1,2)-2*sum(nv==1)*vL(1,2)/dt/L(1,1);
     nv = (v(:,1) < 0).*(r(:,1) < D/2); v(nv==1,1) = -v(nv==1,1); pressure_q(q,1)=pressure_q(q,1)+2*sum(v(nv==1,1))/dt/L(1,2);
     nv = (v(:,2) < 0).*(r(:,2) < D/2); v(nv==1,2) = -v(nv==1,2); pressure_q(q,1)=pressure_q(q,1)+2*sum(v(nv==1,2))/dt/L(1,1);
     
     % Calculate entropy from velocity distribution
    c_vrxy = M^3*floor((v(:,1) - min_vx)/d_vx + 1/2) + M^2*floor((v(:,2) - min_vy)/d_vy + 1/2) + ...
             M*floor((r(:,1) - min_rx)/d_rx) + floor((r(:,2) - min_ry)/d_ry) + 1;
    [n_vrxy(q,1:M^4)] = hist(c_vrxy,(1:M^4));
    
    %calculate temperature
    temperature_q(q,1)= 2*sum(v(:,1).^2+v(:,2).^2)/N;
    
    r = r+v*dt;
    end
    
    volume(t,1)= L(1,1)*L(1,2);
    pressure(t,1)=sum(pressure_q)/Q;
    temperature(t,1)=sum(temperature_q)/Q;
    
    p_vrxy= sum(n_vrxy)/sum(sum(n_vrxy))/d_vrxy;
    entropy(t,1) =sum(-p_vrxy(p_vrxy>0).*log(p_vrxy(p_vrxy>0)))*d_vrxy;
    
    plot(log(volume(1:t,1)),entropy(1:t,1),'r','Linewidth',2);
    plot(log(volume(1:t,1)),temperature(1:t,1),'r','Linewidth',2);
    plot(log(volume(1:t,1)),pressure(1:t,1),'r','Linewidth',2);
    L = L+vL*Q*dt;
    
    f=polyfit(log(volume(1:t,1)),entropy(1:t,1),1);
    [t f]
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    