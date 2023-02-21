// Declare variables
var Lx = 600;
var Ly = 600;
var nx = 128;
var ny = Math.floor(nx * Ly / Lx);
var nghost = 2;
var nnx = nx + 2*nghost;
var nny = ny + 2*nghost;
var nx1 = nghost + 1;
var nx2 = nghost + nx;
var ny1 = nghost + 1;
var ny2 = nghost + ny;
var dt;
var it = 0;
var gam = 5./3.;
var gm1 = gam - 1.0;
var cfl = 0.2;
var ncolors = 255;
var large_number = 1.00e+30;
var nout = 2;
var dmin  = 0.8;
var dmax = 2.2;
var dtnew;
var dfix = 0.8;
var pfix = 1.0e-4;
var icell = 1;
var jcell = 1;

var pClick = 1000.0;
var w0 = 0.1;
var sigma = 0.05/Math.sqrt(2.);

// Arrays

var dens_old = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  dens_old[i] = new Array(nny);
}

var dens_new = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  dens_new[i] = new Array(nny);
}

var momx_old = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  momx_old[i] = new Array(nny);
}

var momx_new = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  momx_new[i] = new Array(nny);
}

var momy_old = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  momy_old[i] = new Array(nny);
}

var momy_new = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  momy_new[i] = new Array(nny);
}

var ener_old = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  ener_old[i] = new Array(nny);
}

var ener_new = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  ener_new[i] = new Array(nny);
}

var update_count = new Array(nnx);
for (var i = 0; i < nnx; i++) {
  update_count[i] = new Array(nny);
}

var red   = new Array(ncolors);
var green = new Array(ncolors);
var blue  = new Array(ncolors);

var varxl = new Array(4);
var varxr = new Array(4);
var varyl = new Array(4);
var varyr = new Array(4);
var fluxx = new Array(4);
var fluxy = new Array(4);

// Gaussian kernel
var kw = 20;
var ngk = 2*kw + 1;
var gkernel = new Array(ngk);
for (var i = 0; i < ngk; i++) {
  gkernel[i] = new Array(ngk);
}
var x = 0.0;
var gsum = 0.0;
var x0 = (ngk+1)*0.5;
var y  = x0*x0/(16.0*Math.log(2.0));
for (var i = 0; i < ngk; i++) {
    for (var j = 0; j < ngk; j++) {
        x = ((i+1)-x0)*((i+1)-x0) + ((j+1)-x0)*((j+1)-x0);
        gkernel[i][j] = Math.exp(-x/y);
        gsum = gsum + gkernel[i][j];
    }
}
for (var i = 0; i < ngk; i++) {
    for (var j = 0; j < ngk; j++) {
       gkernel[i][j] = gkernel[i][j]/gsum;
    }
}


//######################################################################
// Generate canvas
//######################################################################

// Determine cell sizes in pixels
var celldx = Lx / nx;
var celldy = Ly / ny;
var cellid;

document.write("<canvas id=\"myCanvas\" width=\""+Lx+"\" height=\""+Ly+"\" style=\"border:0px;margin:0px;padding:0px;\">Your browser does not support the canvas element.</canvas>");


//######################################################################
// Generate colormap
//######################################################################

for (var i = 0; i < ncolors; i++) {
    var cr = 255.;
    var cg = 255.;
    var cb = 255.;
    if (i < (0.25 * ncolors)) {
       cr = 0;
       cg = Math.floor(4 * i);
    } else if (i < (0.5 * ncolors)) {
       cr = 0;
       cb = Math.floor(1 + 4 * (0.25 * ncolors - i));
    } else if (i < (0.75 * ncolors)) {
       cr = Math.floor(4 * (i - 0.5 * ncolors));
       cb = 0;
    } else {
       cg = Math.floor(1 + 4 * (0.75 * ncolors - i));
       cb = 0;
    }
    red[i]   = cr;
    green[i] = cg;
    blue[i]  = cb;
}























//######################################################################
// Simulation setup
//######################################################################

for (var i = 0; i < nnx; i++) {
    for (var j = 0; j < nny; j++) {
        
        var Y = 1.*(j+0.5)/ny;
        var X = 1.*(i+0.5)/nx;
        var perturb = 0.0;
        if(Math.abs(Y-0.5) < 0.25) { 
          perturb = 1.;
        }

        // Conservative variables
        dens_old[i][j] = (1.0 + perturb);
        momx_old[i][j] = dens_old[i][j] * (-0.5 + perturb);
        momy_old[i][j] = dens_old[i][j] * (w0*Math.sin(4.*Math.PI*X) * ( Math.exp(-Math.pow(Y-0.25,2)/(2*sigma*sigma)) + Math.exp(-Math.pow(Y-0.75,2)/(2*sigma*sigma)) ));
        ener_old[i][j] = 2.5 / gm1 + 0.5*(Math.pow(momx_old[i][j],2) + Math.pow(momy_old[i][j],2)) / dens_old[i][j];
             
        dens_new[i][j] = dens_old[i][j];
        momx_new[i][j] = momx_old[i][j];
        momy_new[i][j] = momy_old[i][j];
        ener_new[i][j] = ener_old[i][j];
        
        update_count[i][j] = 0;
    }
}

// Compute initial timestep
compute_timestep();




// Add click event listener on the canvas
var theCanvas = document.getElementById('myCanvas');
theCanvas.addEventListener("click", getClickPosition, false);

function getClickPosition(e) {
    icell = Math.round(e.clientX/celldx);
    jcell = Math.round(e.clientY/celldy);
    
    // Call energy injection
    inject_energy();
};




// Create animation loop
var mainloop = function() {
    compute_hydro();
    if (it % nout == 0) {
       update_screen();
    }
};

var animFrame = window.requestAnimationFrame       ||
                window.webkitRequestAnimationFrame ||
                window.mozRequestAnimationFrame    ||
                window.oRequestAnimationFrame      ||
                window.msRequestAnimationFrame     ||
                null ;

var recursiveAnim = function() {
    mainloop();
    animFrame( recursiveAnim );
};

// start the mainloop
animFrame( recursiveAnim );


// Inject energy when mouse left button is clicked.
// This is smoothed over a small region using a gaussian kernel
function inject_energy() {
    
    var i1 = Math.max(icell-kw,nx1);
    var i2 = Math.min(icell+kw,nx2);
    var j1 = Math.max(jcell-kw,ny1);
    var j2 = Math.min(jcell+kw,ny2);

    for (var i = i1; i < i2+1; i++) {
        for (var j = j1; j < j2+1; j++) {
            
            ener_old[i][j] = ener_old[i][j] + pClick * gkernel[i-icell+kw][j-jcell+kw] / gm1;
            ener_new[i][j] = ener_old[i][j];
            
            dtcell = dt_cell(i,j);
            dtnew = Math.min(dtnew,dtcell);
            
        }
    }
    
};


//######################################################################
// Classical hydrodynamics routines
//######################################################################

// Compute initial timestep
function compute_timestep() {

    var soundspeed;
    var wsx;
    var wsy;
    var dtx;
    var dty;
    var kinetic_energy;
    var internal_energy;
    
    dt = large_number;
    

    for (var i = nx1-1; i < nx2; i++) {
        for (var j = ny1-1; j < ny2; j++) {
            
            // Compute sound speed inside cell
            kinetic_energy = 0.5 * (momx_new[i][j]*momx_new[i][j] + momy_new[i][j]*momy_new[i][j]) / dens_new[i][j];
            
            internal_energy = ener_new[i][j] - kinetic_energy;
         
            soundspeed = Math.sqrt( gam * gm1 * internal_energy / dens_new[i][j] );
         
            wsx = Math.abs(momx_new[i][j])/dens_new[i][j] + soundspeed;
            wsy = Math.abs(momy_new[i][j])/dens_new[i][j] + soundspeed;
            
            dtx = 1.0 / wsx;
            dty = 1.0 / wsy;
         
            dt = Math.min(dt,dtx,dty);
            
        }
    }
          
  // CFL condition
  dt = cfl * dt
  
};

// Compute timestep in given cell
// Also perform a pressure fix at the same time to keep the number of loops to a minimum
function dt_cell(ii,jj) {
    
    // Pressure fix
    var density = dens_new[ii][jj]
    var kinetic_energy = 0.5 * (momx_new[ii][jj]*momx_new[ii][jj] + momy_new[ii][jj]*momy_new[ii][jj]) / density;
    var internal_energy = ener_new[ii][jj] - kinetic_energy;

    if(density < dfix){
        // maintain the same velocity and pressure in cell whilst resetting density
        dens_new[ii][jj] = dfix
        momx_new[ii][jj] = momx_new[ii][jj] / density * dfix
        momy_new[ii][jj] = momy_new[ii][jj] / density * dfix
        kinetic_energy = 0.5 * (momx_new[ii][jj]*momx_new[ii][jj] + momy_new[ii][jj]*momy_new[ii][jj]) / dfix;
        density = dfix
    }
    
    // internal energy
    if(internal_energy < pfix*kinetic_energy){
        internal_energy = pfix * kinetic_energy
    }

    // update total energy
    ener_new[ii][jj] = internal_energy + kinetic_energy
    
    // Compute sound speed inside cell     
    var soundspeed = Math.sqrt( gam * gm1 * internal_energy / density );
 
    // Wave speeds
    var wsx = Math.abs(momx_new[ii][jj])/density + soundspeed;
    var wsy = Math.abs(momy_new[ii][jj])/density + soundspeed;
    
    var dtx = cfl / wsx;
    var dty = cfl / wsy;
 
    var dtc = Math.min(dtx,dty);
  
    return dtc;
  
};

//######################################################################

// Main hydro loop to compute fluxes and timestep in the entire grid
function compute_hydro() {
    
    it = it + 1; // Increase timestep number
    dtnew = large_number;
    
    // fill ghosts cells
    for (var j = ny1-1; j < ny2; j++) {
        for (var i = 0; i < nghost; i++) {
            dens_old[i][j] = dens_old[nx2-2+i][j];
            momx_old[i][j] = momx_old[nx2-2+i][j];
            momy_old[i][j] = momy_old[nx2-2+i][j];
            ener_old[i][j] = ener_old[nx2-2+i][j];
        }
        for (var i = nx2; i < nnx; i++) {
            dens_old[i][j] = dens_old[i-nx2+3][j];
            momx_old[i][j] = momx_old[i-nx2+3][j];
            momy_old[i][j] = momy_old[i-nx2+3][j];
            ener_old[i][j] = ener_old[i-nx2+3][j];
        }
    }
    
    for (var i = nx1-1; i < nx2; i++) {
        for (var j = 0; j < nghost; j++) {
            dens_old[i][j] = dens_old[i][ny2-2+j];
            momx_old[i][j] = momx_old[i][ny2-2+j];
            momy_old[i][j] = momy_old[i][ny2-2+j];
            ener_old[i][j] = ener_old[i][ny2-2+j];
        }
        for (var j = ny2; j < nny; j++) {
            dens_old[i][j] = dens_old[i][j-ny2+3];
            momx_old[i][j] = momx_old[i][j-ny2+3];
            momy_old[i][j] = momy_old[i][j-ny2+3];
            ener_old[i][j] = ener_old[i][j-ny2+3];
        }
    }
    
    
    for (var i = nx1-1; i < nx2+1; i++) {
        for (var j = ny1-1; j < ny2+1; j++) {
         
            // Compute values at the cell interfaces using gradients
            compute_interface_values([dens_old[i-2][j],dens_old[i-1][j],dens_old[i][j],dens_old[i+1][j]],[momx_old[i-2][j],momx_old[i-1][j],momx_old[i][j],momx_old[i+1][j]],[momy_old[i-2][j],momy_old[i-1][j],momy_old[i][j],momy_old[i+1][j]],[ener_old[i-2][j],ener_old[i-1][j],ener_old[i][j],ener_old[i+1][j]],varxl,varxr)
            
            compute_interface_values([dens_old[i][j-2],dens_old[i][j-1],dens_old[i][j],dens_old[i][j+1]],[momy_old[i][j-2],momy_old[i][j-1],momy_old[i][j],momy_old[i][j+1]],[momx_old[i][j-2],momx_old[i][j-1],momx_old[i][j],momx_old[i][j+1]],[ener_old[i][j-2],ener_old[i][j-1],ener_old[i][j],ener_old[i][j+1]],varyl,varyr)
            
            // Compute fluxes using Riemann solver
            riemann_solver(varxl,varxr,fluxx);
            riemann_solver(varyl,varyr,fluxy);
            
            // Update conservative variables on each side of the interface
            dens_new[i  ][j] = dens_new[i  ][j] + fluxx[0]*dt;
            dens_new[i-1][j] = dens_new[i-1][j] - fluxx[0]*dt;
            momx_new[i  ][j] = momx_new[i  ][j] + fluxx[1]*dt;
            momx_new[i-1][j] = momx_new[i-1][j] - fluxx[1]*dt;
            momy_new[i  ][j] = momy_new[i  ][j] + fluxx[2]*dt;
            momy_new[i-1][j] = momy_new[i-1][j] - fluxx[2]*dt;
            ener_new[i  ][j] = ener_new[i  ][j] + fluxx[3]*dt;
            ener_new[i-1][j] = ener_new[i-1][j] - fluxx[3]*dt;
            
            dens_new[i][j  ] = dens_new[i][j  ] + fluxy[0]*dt;
            dens_new[i][j-1] = dens_new[i][j-1] - fluxy[0]*dt;
            momy_new[i][j  ] = momy_new[i][j  ] + fluxy[1]*dt;
            momy_new[i][j-1] = momy_new[i][j-1] - fluxy[1]*dt;
            momx_new[i][j  ] = momx_new[i][j  ] + fluxy[2]*dt;
            momx_new[i][j-1] = momx_new[i][j-1] - fluxy[2]*dt;
            ener_new[i][j  ] = ener_new[i][j  ] + fluxy[3]*dt;
            ener_new[i][j-1] = ener_new[i][j-1] - fluxy[3]*dt;
            
            // Count how many faces have been updated
            update_count[i  ][j] = update_count[i  ][j] + 2;
            update_count[i-1][j] = update_count[i-1][j] + 1;
            update_count[i][j-1] = update_count[i][j-1] + 1;
            
            // If all faces have been updates, the cell will not be revisited, so we can
            // finalise the update and compute the new timestep.
            // we only need to check in i-1, not j-1 because of the order of the i,j loops
            if(update_count[i-1][j] == 4){
                
                dtcell = dt_cell(i-1,j);
                dtnew = Math.min(dtnew,dtcell);
                
                dens_old[i-1][j] = dens_new[i-1][j];
                momx_old[i-1][j] = momx_new[i-1][j];
                momy_old[i-1][j] = momy_new[i-1][j];
                ener_old[i-1][j] = ener_new[i-1][j];
                
                update_count[i-1][j] = 0;
                
            }
            
        }
    }
    
    dt = dtnew;
    
};
     
//######################################################################

function riemann_solver(ucl,ucr,flux) {

    var dl,pl,ul,vl,el;
    var dr,pr,ur,vr,er;
    var cfastl,cfastr;
    var Cmax;
    
    dl = ucl[0]
    ul = ucl[1] / ucl[0]
    vl = ucl[2] / ucl[0]
    el = ucl[3]
    pl = ( ucl[3] - 0.5*(ucl[1]*ucl[1]+ucl[2]*ucl[2])/ucl[0] ) * gm1
    
    dr = ucr[0]
    ur = ucr[1] / ucr[0]
    vr = ucr[2] / ucr[0]
    er = ucr[3]
    pr = ( ucr[3] - 0.5*(ucr[1]*ucr[1]+ucr[2]*ucr[2])/ucr[0] ) * gm1
    
    // Find the largest eigenvalues in the normal direction to the interface
    cfastl = Math.sqrt(gam*pl/dl)
    cfastr = Math.sqrt(gam*pr/dr)
    
    // Compute max wave speed
    Cmax = Math.max(Math.abs(ul),Math.abs(ur)) + Math.max(Math.abs(vl),Math.abs(vr)) + Math.max(cfastl,cfastr)
    
    flux[0] = 0.5*(dl*ul + dr*ur) + 0.5*Cmax*(dl - dr);
    flux[1] = 0.5*((dl*ul*ul + pl) + (dr*ur*ur + pr)) + 0.5*Cmax*((dl*ul) - (dr*ur));
    flux[2] = 0.5*((dl*ul*vl) + (dr*ur*vr)) + 0.5*Cmax*((dl*vl) - (dr*vr));
    flux[3] = 0.5*(((el + pl) * ul) + ((er + pr) * ur)) + 0.5*Cmax*(el - er);

};



// Compute interface values with slopes
function compute_interface_values(dens,momx,momy,ener,varl,varr) {

    var dlft;
    var drgt;
    var dcen;
    var dsgn;
    var dlim;
    var av;
    
    // Left state
    
    dlft = 2.0*(dens[1]-dens[0]);
    drgt = 2.0*(dens[2]-dens[1]);
    dcen = 0.5*(dens[2]-dens[0]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varl[0] = dens[1] + av;
    
    dlft = 2.0*(momx[1]-momx[0]);
    drgt = 2.0*(momx[2]-momx[1]);
    dcen = 0.5*(momx[2]-momx[0]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varl[1] = momx[1] + av;
    
    dlft = 2.0*(momy[1]-momy[0]);
    drgt = 2.0*(momy[2]-momy[1]);
    dcen = 0.5*(momy[2]-momy[0]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varl[2] = momy[1] + av;
    
    dlft = 2.0*(ener[1]-ener[0]);
    drgt = 2.0*(ener[2]-ener[1]);
    dcen = 0.5*(ener[2]-ener[0]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varl[3] = ener[1] + av;
    
    
    // Right state
    
    dlft = 2.0*(dens[2]-dens[1]);
    drgt = 2.0*(dens[3]-dens[2]);
    dcen = 0.5*(dens[3]-dens[1]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varr[0] = dens[2] - av;
    
    dlft = 2.0*(momx[2]-momx[1]);
    drgt = 2.0*(momx[3]-momx[2]);
    dcen = 0.5*(momx[3]-momx[1]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varr[1] = momx[2] - av;
    
    dlft = 2.0*(momy[2]-momy[1]);
    drgt = 2.0*(momy[3]-momy[2]);
    dcen = 0.5*(momy[3]-momy[1]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varr[2] = momy[2] - av;
    
    dlft = 2.0*(ener[2]-ener[1]);
    drgt = 2.0*(ener[3]-ener[2]);
    dcen = 0.5*(ener[3]-ener[1]);
    dsgn = Math.sign(dcen);
    dlim = Math.min(Math.abs(dlft),Math.abs(drgt));
    if((dlft*drgt) <= 0.0){
        dlim = 0.0;
    }
    av = 0.5*dsgn*Math.min(dlim,Math.abs(dcen));
    varr[3] = ener[2] - av;
    
    
    
};




//######################################################################
// Update graphics
//######################################################################

function update_screen() {

    var canvas = document.getElementById("myCanvas");
    var ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height); // Clear the canvas and start fresh

    var icol;
  
    for (var i = nx1-1; i < nx2; i++) {
        for (var j = ny1-1; j < ny2; j++) {
            
            icol = Math.min(Math.max(Math.floor(((dens_new[i][j]-dmin)/(dmax-dmin))*(ncolors-1)),0),ncolors-1)
            theString = "rgb("+String(red[icol])+","+String(green[icol])+","+String(blue[icol])+")";
            ctx.fillStyle = theString;
            ctx.fillRect(i*celldx,j*celldy,celldx,celldy);
        }
    }
    
};
