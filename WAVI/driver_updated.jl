

using WAVI


function coupled_driver()

### Grid and boundary conditions ###
nx = 468
ny  = 398
nσ = 12
#swapped around chekc still right.
x0 = -838500.0 
y0 = -1792500.0
dx = 2000.0
dy = 2000.0

Cxl = 44
Cxu = 320
Cyl = 190
Cyu = 375

x_w=0  
x_e=0     
y_s=0
y_n=0
    if Cxl > 1 || Cxu < nx || Cyl > 1 || Cyu < ny  
     if Cxl > 1 
#      x_w=x_w+1
     end  
     if Cxu <  nx
      x_e=x_e+1
     end     
     if Cyl > 1 
      y_s=y_s+1
     end   
     if Cyu <  ny
#      y_n=y_n+1
     end
    end


#timestepping parameters
niter0 = 0 #CHANGE ME FOR A PICKUP (!! must be niter0 !!)
#step_thickness = false
step_thickness = true  #overall stepping thickness
step_haf = true       #in Ocean /ice area
dt_s = unset
dt = dt_s/(3600*24*365)
end_time_s = unset
end_time =end_time_s/(3600*24*365)
chkpt_freq_s = unset
chkpt_freq = chkpt_freq_s/(3600*24*365)
pchkpt_freq_s = unset
pchkpt_freq = pchkpt_freq_s/(3600*24*365)
Couple_run_end = unset
dt_coup = unset
dt_coup = dt_coup/1


if niter0 == 0

h_mask=Array{Float64}(undef, nx, ny);
read!("Inverse_2km_h_mask_clip_rot_02_BedmachineV3_FULL_stripe_fix.bin", h_mask)
#read!("Inverse_2km_h_mask_clip_rot_BedmachineV2.bin", h_mask)
h_mask.=ntoh.(h_mask)
#swapped as roatted
u_iszero=Array{Float64}(undef,nx+1,ny);
read!("Inverse_2km_viszero_clip_rot_02_BedmachineV3_FULL_stripe_fix.bin",u_iszero)
#read!("Inverse_2km_viszero_clip_rot_BedmachineV2.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("Inverse_2km_uiszero_clip_rot_02_BedmachineV3_FULL_stripe_fix.bin",v_iszero)
#read!("Inverse_2km_uiszero_clip_rot_BedmachineV2.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("Inverse_2km_sigma_grid_BedmachineV3_FULL_stripe_fix.bin",sigma_grid)
#read!("Inverse_2km_sigma_grid_BedmachineV2.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)


grid = Grid(nx = nx, 
            ny = ny,   
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            h_mask = h_mask, 
            u_iszero = u_iszero, 
            v_iszero = v_iszero,
            σ = sigma_grid,
            Cxl = Cxl,
            Cxu = Cxu,
            Cyl = Cyl,
            Cyu = Cyu)

#Bed 
bed_elevation= Array{Float64}(undef, nx, ny);
read!("Inverse_2km_bed_clip_noNan_rot_BedmachineV3_FULL_edit_all_280_taper_6000_min_140.bin", bed_elevation)
bed_elevation .= ntoh.(bed_elevation)

#solver parameters
maxiter_picard = 60
tol_picard = 1e-3
solver_params = SolverParams(maxiter_picard = maxiter_picard,
                        tol_picard = tol_picard)

#initial conditions
starting_thickness= Array{Float64}(undef, nx, ny);
read!("Inverse_2km_thickness_clip_noNan_02_rot_BedmachineV3_FULL_stripe_fix_rel_4000_2.bin", starting_thickness)
starting_thickness .= ntoh.(starting_thickness)


viscosity=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_viscosity3D_clip_noNan_rot_02_BedmachineV3_FULL_stripe_fix.bin",viscosity)
viscosity.=ntoh.(viscosity)

# swap as rotated

initial_guess_u=Array{Float64}(undef, nx + 1, ny);
read!("Inverse_2km_v_velocs_clip_noNan_rot_02_BedmachineV3_FULL_stripe_fix_rel_4000.bin", initial_guess_u)
initial_guess_u .= ntoh.(initial_guess_u)


initial_guess_v=Array{Float64}(undef, nx, ny + 1);
read!("Inverse_2km_u_velocs_clip_noNan_rot_02_BedmachineV3_FULL_stripe_fix_rel_4000.bin", initial_guess_v)
initial_guess_v .= ntoh.(initial_guess_v)

#

temp=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_3Dtemp_clip_noNan_rot_02_BedmachineV3_FULL_stripe_fix.bin",temp)
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_damage3D_clip_noNan_rot_02_BedmachineV3_FULL_stripe_fix.bin",damage)
damage.=ntoh.(damage)

#
 
weertman_c=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_WeertmanC_clip_adjusted_noNan_rot_02_BedmachineV3_FULL_stripe_fix.bin",weertman_c)
weertman_c.=ntoh.(weertman_c)

accumulation_rate=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_accumulation_clip_noNan_rot_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

initial_conditions = InitialConditions(initial_thickness = starting_thickness,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage,
                                        initial_u_veloc = initial_guess_u,
                                        initial_v_veloc = initial_guess_v)

params = Params(accumulation_rate = accumulation_rate,
                weertman_c = weertman_c,
		step_haf = step_haf)

#make the model
model = Model(grid = grid,
                    bed_elevation = bed_elevation, 
                    params = params, 
                    solver_params = solver_params,
                    initial_conditions = initial_conditions)

end

if niter0 == 0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                        dt = dt, 
                                        end_time = end_time, 
                                        chkpt_freq = chkpt_freq, 
                                        pchkpt_freq = pchkpt_freq,
                                        step_thickness = step_thickness)
else
timestepping_params_new = TimesteppingParams(niter0 = niter0,
                                        dt = dt,
                                        end_time = end_time,
                                        chkpt_freq = chkpt_freq,
                                        pchkpt_freq = pchkpt_freq,
                                        step_thickness = step_thickness)

end

if niter0 == 0

#### output parameters ###

outputs = (h = model.fields.gh.h,
           uh = model.fields.gh.u,
           dhdt = model.fields.gh.dhdt,
           u = model.fields.gu.u,
           vh = model.fields.gh.v,
	   dep_vis= model.fields.gh.ηav,
           W_c= model.fields.gh.weertman_c,
           av_sp= model.fields.gh.av_speed,
           beta2= model.fields.gh.β,
           glen_b= model.fields.g3d.glen_b,
           damage= model.fields.g3d.Φ,
           v = model.fields.gv.v,
           us = model.fields.gh.us,
           vs = model.fields.gh.vs,
           ub = model.fields.gh.ub,
           vb = model.fields.gh.vb,
           haf = model.fields.gh.haf,
           gh_frac = model.fields.gh.grounded_fraction,
           gu_frac = model.fields.gu.grounded_fraction,
           gv_frac = model.fields.gv.grounded_fraction,
           h_mask = model.fields.gh.mask)


output_freq_s = unset
output_freq = output_freq_s/(3600*24*365)
output_params = OutputParams(outputs = outputs,
                             output_freq = output_freq,
                             output_format = "mat",
#                             zip_format = "nc",
                             dump_vel = true,
                             PC_east = true,
                             PC_south = true,
			     dt_coup = dt_coup,
			     Output_vel = true,
                             Output_Hb = true,
                             Output_float = true)

end

if niter0 > 0

n_iter_string =  lpad(timestepping_params_new.niter0, 10, "0"); #filename as a string with 10 digits
        println("detected niter0 > 0 (niter0 = $(timestepping_params_new.niter0)). Looking for pickup...")
        simulation=load(string("PChkpt_",n_iter_string, ".jld2"), "simulation")

simulation = @set simulation.timestepping_params = timestepping_params_new

else

### set up the simulation ###
#!! This does a pickup if niter0 >0, does not rebuild the model unless flag activated !!
simulation = Simulation(model = model,
                        timestepping_params = timestepping_params,
                        output_params = output_params)

end






### set the thickness and float fields if we're doing a pickup
if niter0 > 0
#set thickness

    thickness_file= Array{Float64}(undef, (simulation.model.grid.Cxu - simulation.model.grid.Cxl +1 + 2 + x_e +x_w), (simulation.model.grid.Cyu- simulation.model.grid.Cyl +1 +2 + y_s +y_n));
    clock_time_C=Int(round((simulation.clock.time*(3600*24*365))/dt_coup))
    Thickness_file_string = string("Streamice_Thickness_out", lpad(clock_time_C, 10,"0") ,".data")
    read!(Thickness_file_string, thickness_file)
#    read!(string("Streamice_Thickness_out.data"), thickness_file)
    thickness_file_noBc=thickness_file[2+x_w:end-1-x_e,2+y_s:end-1-y_n];
    thickness_file_noBc .= ntoh.(thickness_file_noBc);


    h_temp=deepcopy(simulation.model.fields.gh.h);
    h_temp[simulation.model.grid.Cxl:simulation.model.grid.Cxu,simulation.model.grid.Cyl:simulation.model.grid.Cyu]=thickness_file_noBc;
    simulation.model.fields.gh.h[simulation.model.fields.gh.mask]=h_temp[simulation.model.fields.gh.mask];



end

#set haf fields 
if step_haf == false

#    haf_file= Array{Float64}(undef, (simulation.model.grid.Cxu - simulation.model.grid.Cxl +1 + 2 + x_e +x_w), (simulation.model.grid.Cyu- simulation.model.grid.Cyl +1 +2 + y_s +y_n));
#    clock_time_C=Int(round((simulation.clock.time*(3600*24*365))/dt_coup))
#    haf_file_string = string("Shelfice_haf_out", lpad(clock_time_C, 10,"0") ,".data")
#    read!(haf_file_string, haf_file)
#    haf_file_noBc=haf_file[2+x_w:end-1-x_e,2+y_s:end-1-y_n];
#    haf_file_noBc .= ntoh.(haf_file_noBc);
#    f_temp=simulation.model.fields.gh.haf;
#    f_temp[simulation.model.grid.Cxl:simulation.model.grid.Cxu,simulation.model.grid.Cyl:simulation.model.grid.Cyu]=haf_file_noBc

#    f_temp_oc=f_temp[simulation.model.grid.Cxl:simulation.model.grid.Cxu,simulation.model.grid.Cyl:simulation.model.grid.Cyu];
#    h_temp=simulation.model.fields.gh.mask;
#    h_temp_oc=h_temp[simulation.model.grid.Cxl:simulation.model.grid.Cxu,simulation.model.grid.Cyl:simulation.model.grid.Cyu];
#    f_temp_oc[h_temp_oc]=haf_file_noBc[h_temp_oc];
#    f_temp[simulation.model.grid.Cxl:simulation.model.grid.Cxu,simulation.model.grid.Cyl:simulation.model.grid.Cyu]=f_temp_oc;
    
end    



### run the simulation ###
run_simulation!(simulation)


if end_time_s == Couple_run_end

  nc_name_full = string(simulation.output_params.output_path, simulation.output_params.prefix, ".nc")
  WAVI.make_ncfile(simulation.output_params.output_format, simulation.output_params.output_path, nc_name_full, simulation.output_params.prefix)
println("WAVI output zipped")

end


return nothing

end

coupled_driver()
