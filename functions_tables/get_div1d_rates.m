function proc = div1d_post_processing(output)

    amu = 1.66e-27;
    e = 1.602e-19;
    
    m_i = 2; %ions
    m_a = 2; %atoms
    m_m = 4; %molecules

    % get data from div1d input/output
    v = output.velocity(end,:);
    ne = output.density(end,:);
    Te = output.temperature(end,:);
    na = output.neutral_density(end,:);
    Ti = Te; 

    try 
        vn = output.neutral_velocity(end,:);
    catch
        vn = zeros(size(v));
    end
    try
        nm = output.molecule_density(end,:);
    catch
        nm = zeros(size(v));
    end
        
    % reduced masses: 
    mu_a = amu*(m_a*m_i/(m_a+m_i));
    mu_m = amu*(m_a*m_i/(m_a+m_i));
    
    % get relevant rate coefficients
    rate_cx_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.2', '3.1.8');
    rate_E_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.3', '0.1D');
    rate_EIR_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.1.8');
    rate_ion_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.1.5');
    
    % Get relevant molecular rate coefficients
    rate_cx_mol_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.2', '3.2.3');
    rate_ion_mol_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.2.9');
    rate_diss_mol_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.2.5');
    rate_E_mol_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.3', '0.3D');
    rate_da_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.2', '2.2.17');
    rate_ene_m_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.10', '2.2.h2c');
    
    % Get relevant H2+ rate coefficients
    rate_diss_rec_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.2.14');
    rate_diss_H2plus_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.2.12');
    rate_diss_ion_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '2.2.11');
    
    % Get relevant H- rate coefficients
    rate_cx_Hmin_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '7.2.3a');
    rate_ion_Hmin_coeffs = get_coeffs('Data/rates/amjuel.tex', 'H.4', '7.2.3b');
     
    % calculate losses along divertor leg
    cx = zeros(length(v),1);
    cx_mom_loss = zeros(length(v),1);
    el_mom_loss = zeros(length(v),1);
    el_mol_mom_loss = zeros(length(v),1);
    MAR = zeros(length(v),1);
    MAI = zeros(length(v),1);
    MAD = zeros(length(v),1);

    ion = zeros(length(v),1);
    EIR = zeros(length(v),1);
    MAR_H2plus = zeros(length(v),1);
    MAD_H2plus = zeros(length(v),1);
    MAR_Hmin = zeros(length(v), 1);
    MAD_Hmin = zeros(length(v), 1);
    MAR_mom_loss = zeros(length(v),1);

    ion_m =  zeros(length(v), 1);
    diss = zeros(length(v), 1);
    cx_m = zeros(length(v), 1);
    da_m = zeros(length(v), 1);
    ene_m = zeros(length(v), 1);
    
    Qd = zeros(length(v), 1);
    Qdi = zeros(length(v), 1);
    Qcx = zeros(length(v), 1);
    Qi = zeros(length(v), 1);
    Qdr = zeros(length(v), 1);
    
    
    
    for i = 1:length(v)
        
        % 'Beam' Energy
        E_a = amu*m_a/(2*e)*(v(i)-vn(i))^2;
        E_m = amu*m_m/(2*e)*v(i)^2;
    
        
        % Momentum loss through charge exchange
        cx(i) = ne(i)*na(i)*eval_1D(rate_cx_coeffs, Ti(i)/m_i);
        cx_mom_loss(i) = amu*m_i*(v(i)-vn(i))*cx(i);
        
        
        % Calculate molecular reaction rates
        rate_cx_mol = eval_1D(rate_cx_mol_coeffs, max(Ti(i)*1.6726e-27/3.3436e-27, 0.1));
        rate_ion_mol = eval_2D_Tne(rate_ion_mol_coeffs, Te(i), ne(i));
        rate_diss_mol = eval_2D_Tne(rate_diss_mol_coeffs, Te(i), ne(i));
        rate_da_mol = eval_1D(rate_da_coeffs, Te(i));
        rate_ene_mol = eval_2D_Tne(rate_ene_m_coeffs, max(Te(i), 0.1), ne(i));
    
        % Calculate H2plus rates
        rate_diss_rec = eval_2D_Tne(rate_diss_rec_coeffs, Te(i), ne(i));
        rate_diss_H2plus = eval_2D_Tne(rate_diss_H2plus_coeffs, Te(i), ne(i));
        rate_diss_ion = eval_2D_Tne(rate_diss_ion_coeffs, Te(i), ne(i));
    
        % Calculate H- rates
        rate_cx_Hmin = eval_2D_Tne(rate_cx_Hmin_coeffs, Ti(i)/m_a, ne(i));
        rate_ion_Hmin = eval_2D_Tne(rate_ion_Hmin_coeffs, Ti(i)/m_a, ne(i));
        
        f_cx = rate_cx_mol / (rate_cx_mol + rate_ion_mol); % fraction of cx H2+ generation
        
        MAR_H2plus_r = (rate_cx_mol+rate_ion_mol)/(rate_diss_rec+rate_diss_H2plus+rate_diss_ion)*rate_diss_rec*f_cx;
        MAD_H2plus_r = (rate_cx_mol+rate_ion_mol)/(rate_diss_rec+rate_diss_H2plus+rate_diss_ion)*(f_cx*rate_diss_mol + (1-f_cx)*rate_diss_rec); 
        MAI_H2plus_r = (rate_cx_mol+rate_ion_mol)/(rate_diss_rec+rate_diss_H2plus+rate_diss_ion)*(f_cx*rate_diss_ion + (1-f_cx)*(rate_diss_mol+2*rate_diss_ion)); 
        MAR_Hmin_r = rate_da_mol/(rate_cx_Hmin+rate_ion_Hmin)*rate_cx_Hmin; 
        MAD_Hmin_r = rate_da_mol/(rate_cx_Hmin+rate_ion_Hmin)*rate_ion_Hmin;
    
        nHplus = nm(i)*(rate_cx_mol+rate_ion_mol)/(rate_diss_rec+rate_diss_H2plus+rate_diss_ion);
        nHmin = nm(i)*rate_da_mol/(rate_cx_Hmin+rate_ion_Hmin);
    
        Qd(i) = ne(i)*nHplus*2.69*rate_diss_H2plus;
        Qdi(i) = nHplus*ne(i)*rate_diss_ion * (2.69+13.6);
	    Qcx(i) = nHmin*ne(i)*rate_cx_Hmin*0.754;
	    Qi(i) = nHmin*ne(i)*rate_ion_Hmin*(0.754+13.6);
        Qdr(i) = nHplus*ne(i)*rate_diss_rec * (2.69-13.6);
    
        MAR_H2plus(i) = ne(i)*nm(i)*MAR_H2plus_r;
        MAD_H2plus(i) = ne(i)*nm(i)*MAD_H2plus_r;
        MAI(i) = ne(i)*nm(i)*MAI_H2plus_r;

        MAR_Hmin(i) = ne(i)*nm(i)*MAR_Hmin_r;
        MAD_Hmin(i) = ne(i)*nm(i)*MAD_Hmin_r;
        
        MAR(i) = MAR_H2plus(i)+MAR_Hmin(i);
        MAD(i) = MAD_H2plus(i)+MAD_Hmin(i);

        ion_m(i) =  ne(i)*nm(i)*rate_ion_mol;
        diss(i) = ne(i)*nm(i)*rate_diss_mol;
        cx_m(i) = ne(i)*nm(i)*rate_cx_mol;
        da_m(i) = ne(i)*nm(i)*rate_da_mol;
        ene_m(i) = ne(i)*nm(i)*rate_ene_mol; 
    
        EIR(i) = ne(i)^2*eval_2D_Tne(rate_EIR_coeffs,max(Te(i),0.1), ne(i));
        
        if E_a>0.1
            el_mom_loss(i) = mu_a * (v(i)-vn(i))*ne(i)*na(i)*eval_2D_TE(rate_E_coeffs, Ti(i)/m_i, E_a);
        else
            el_mom_loss(i) = 0; 
        end
    
        if E_m>0.1
            el_mol_mom_loss(i) = mu_m*v(i)*ne(i)*nm(i)*eval_2D_TE(rate_E_mol_coeffs, Ti(i)/m_i, E_m);
        else
            el_mol_mom_loss(i) = 0;
        end
        
    %     el_mom_loss(i) = ne(i)*na(i)*eval_2D_TE(rate_E_coeffs, max(Ti(i)/m_i, 0.1), max(E_a, 0.1));
    %     el_mol_mom_loss(i) = ne(i)*nm(i)*eval_2D_TE(rate_E_mol_coeffs, max(Ti(i)/m_i, 0.1), max(E_m, 0.1));
        
        % Ion sources
        rate_ion = eval_2D_Tne(rate_ion_coeffs, max(Te(i), 0.1), ne(i));
        ion(i) = ne(i)*na(i)*rate_ion;
    
        MAR_mom_loss(i) = amu*m_i*v(i)*ne(i)*nm(i)*MAR(i);
    
    end
    particles.MAR.value = MAR;
    particles.MAD.value = MAD;
    particles.MAI.value = MAI;
    particles.MAR_H2plus.value = MAR_H2plus;
    particles.MAD_H2plus.value = MAD_H2plus;
    particles.MAR_Hmin.value = MAR_Hmin;
    particles.MAD_Hmin.value = MAD_Hmin;
    particles.EIR.value = EIR; 
    particles.ion_m.value = ion_m;
    particles.diss.value = diss;
    particles.cx_m.value = cx_m;
    particles.da_m.value = da_m;
    particles.ion.value = ion; 

    energy.ene_m.value = ene_m;
    energy.Qd.value = Qd; 
    energy.Qdr.value = Qdr; 
    energy.Qdi.value = Qdi; 
    energy.Qi.value = Qi;
    energy.Qcx.value = Qcx; 

    momentum.MAR.value = MAR_mom_loss;
    momentum.cx.value = cx_mom_loss;
    momentum.EC_n.value = el_mom_loss;
    momentum.EC_m.value = el_mol_mom_loss;

    particles.MAR.label = 'Molecular activated recombination';
    particles.MAD.label = 'Molecular activated dissociation';
    particles.MAI.label = 'Molecular activated ionization';
    particles.MAR_H2plus.label = 'MAR through D2+';
    particles.MAD_H2plus.label = 'MAD through D2+';
    particles.MAR_Hmin.label = 'MAR through D-';
    particles.MAD_Hmin.label = 'MAD through D-';
    particles.EIR.label = 'Electron ion recombination'; 
    particles.ion_m.label = 'Ionization of molecules';
    particles.diss.label = 'Dissociation';
    particles.cx_m.label = 'Molecular charge exchange';
    particles.da_m.label = 'Dissociative attachment';
    particles.ion.label = 'Electron impact ionization of atoms';  

    energy.ene_m.label = 'Energy loss due to dissociation/ionization';
    energy.Qd.label = 'Dissociation D2+'; 
    energy.Qdr.label = 'Dissociative recombination'; 
    energy.Qdi.label = 'Dissociative ionization D2+'; 
    energy.Qi.label = 'Ionization D-';
    energy.Qcx.label = 'Charge exchange D-'; 

    momentum.MAR.label = 'Molecular activated recombination';
    momentum.cx.label = 'Atomic charge exchange';
    momentum.EC_n.label = 'Elastic ion-atom collisions';
    momentum.EC_m.label = 'Elastic ion-molecule collisions';

    fn = fieldnames(particles);
    for i=1:length(fn)
        fn1 = fn{i};
        particles.(fn1).unit = 'm^{-3}s^{-1}'; 
    end

    fn = fieldnames(momentum);
    for i=1:length(fn)
        fn1 = fn{i};
        momentum.(fn1).unit = 'kg m^{-2}s^{-1}'; 
    end

    fn = fieldnames(energy);
    for i=1:length(fn)
        fn1 = fn{i};
        energy.(fn1).unit = 'eV m^{-3}s^{-1}'; 
    end
    proc.particles = particles;
    proc.energy = energy; 
    proc.momentum = momentum; 

end
