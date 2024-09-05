
//---------------------------------//
//     Particle Life simulation    //
//     v0.21a  04/08/2024          //
//                                 //
//                                 //
//     by cienzorama@gmail.com     //
//---------------------------------//

n=1500;
pt=6;

matcol=[];
matcola=[];
matR1=[];
matR2=[];

pcol=[];
px=[];vx=[];ax=[];
py=[];vy=[];ay=[];

K=0.05;
fr=0.85;
dt=1;
dtd2=0.5*dt;

wh=0;wv=0;
cnt=0;

mat_list=[];
mat_n_list=[];
mat_n_list_zeros=[];

mx=0;my=0;
mx_wh=0;my_wv=0;
global_Rmax=250;
global_scale=0.5;

start_time=0;
global_iter=0;
global_time=0;
maxR2dR1=-1;

dc=10;

//Función general de configuración
function setup() {

  //Tamaño del lienzo
  wh=windowWidth-10;
  wv=windowHeight-10;
  whd2=wh*0.5;
  wvd2=wv*0.5;
  whd2n=-wh*0.5;
  wvd2n=-wv*0.5;

  //Tamaños de las matrices de celdas
  mx=Math.round(wh/(global_Rmax*global_scale));
  my=Math.round(wv/(global_Rmax*global_scale));

  mx_wh=mx/wh;
  my_wv=my/wv;

  createCanvas(wh,wv);
  angleMode(DEGREES);

  let force_map_text = createP("<span style='color:white'> Force Map </span>");
  force_map_text.position((pt+4)*dc, (pt+2)*dc/2);
  let R1R2_map_text = createP("<span style='color:white'> R2/R1 Map </span>");
  R1R2_map_text.position((pt+4)*dc, 1.5*(pt+2)*dc+dc);

  //Generación inicial de las 
  //propiedades de las partículas y
  //listas de partículas
  init();
  init_lista();

  start_time = new Date().getTime();
  
}

//Bucle general de representación
function draw() {
  background(0, 0, 0);

  cnt++;
  //Reinicia las reglas de interacción cada 300 frames y
  //calcula el tiempo transcurrido entre cada reinicio
  if (cnt==300) {

    global_iter++;
    var end_time = new Date().getTime();
    var delta_time = end_time - start_time;
    global_time+=delta_time;
    console.log("Iteracion: ",global_iter,
                "Current Time: ",delta_time,
                "Mean Time: ", global_time/global_iter);
    start_time=end_time*1;

    matriz_de_interaccion();
    cnt=0;
  }

  //Integración numérica en cada paso
  integraLF();

  //Representación de las partículas
  noStroke();
  for (var i=0;i<n;i++) {

    fill(Scolors[pcol[i]]);
    circle(px[i], py[i], 5);

  }

  //Representación de fuerzas y radios
  stroke(255);
  var delta_pos=(pt+3)*dc;
  
  var xc=0;var yc=dc;
  for (var i=-2;i<pt;i++) {

    xc+=dc;
    yc=dc;

    for (var j=-2;j<pt;j++) {

      yc+=dc; 

      if (j==-2 && (i==-2 || i==-1)) continue;
      if (j==-2) {
        fill(Scolors[i]);
        circle(xc+dc/2,yc+dc/2,dc*0.75);
        circle(xc+dc/2,yc+dc/2+delta_pos,dc*0.75);
        continue;
      }
      if (j==-1) continue;

      if (i==-2) {
        fill(Scolors[j]);
        circle(xc+dc/2,yc+dc/2,dc*0.75);
        circle(xc+dc/2,yc+dc/2+delta_pos,dc*0.75);
        continue;
      }
      if (i==-1) continue;

      var val1=Math.floor(map(matcol[i][j],-1,1,0,255));
      fill(jet_colormap[val1]);
      rect(xc,yc,dc,dc);

      var val2=Math.floor(map(matR2[i][j]/matR1[i][j],0,maxR2dR1,0,255));
      fill(jet_colormap[val2]);
      rect(xc,yc+delta_pos,dc,dc);
      
    }
  }
  //

}

//Integrador numérico Leap-Frog (E~O2)
function integraLF() {

  for (var i=0;i<n;i++) {

    vx[i]+=dtd2*ax[i];
    vy[i]+=dtd2*ay[i];

    //Actualización de la posición y aplicación
    //de las condiciones de contorno periódicas
    px[i]=(px[i]+dt*vx[i]+wh) % wh;
    py[i]=(py[i]+dt*vy[i]+wv) % wv;

  }

  //Cálculo de la fuerza
  //de interacción
  fuerza2();

  for (var i=0;i<n;i++) {

    vx[i]+=dtd2*ax[i];
    vy[i]+=dtd2*ay[i];

    //Fricción ficticia
    vx[i]*=fr;
    vy[i]*=fr;

  }

}

//Evaluación de la fuerza de
//interacción entre partículas
function fuerza2() {

  var ffid=0;
  var ff=0;
  var xi=0;var yi=0;
  var dtx=0;var dty=0;
  var i2=0;var j2=0;
  var dx=0;var dy=0; var d=0;
  var R1=0;var R2=0; var id=0;
  var Rm=0; var R1dR2=0;

  //Borra el número de elementos de cada celda
  // de la lista de partículas
  mat_n_list=structuredClone(mat_n_list_zeros);

  //Actualiza la lista de elementos existentes
  // en cada celda
  for (var i=0;i<n;i++) {

    //Obtiene la celda correspondiente
    //a la partícula i-ésima
    xi=Math.floor(px[i]*mx_wh);
    yi=Math.floor(py[i]*my_wv);

    //Añade la partícula a la celda correspondiente
    mat_list[xi][yi][ mat_n_list[xi][yi] ] = i;
    (mat_n_list[xi][yi])++;

  }
  
  for (var i=0;i<n;i++) {

    //Borra las interacciones
    //previamente calculadas
    ax[i]=0;ay[i]=0;

    //Obtiene la celda correspondiente
    //a la partícula i-ésima
    xi=Math.floor(px[i]*mx_wh);
    yi=Math.floor(py[i]*my_wv);

    //Busca partículas en la celda actual
    // y en las 8 inmediatamente adyascetes
    for (var i3=xi-1;i3<xi+2;i3++) {
      for (var j3=yi-1;j3<yi+2;j3++) {

        dtx=0;dty=0;
        i2=i3;j2=j3;
        
        //Condiciones de contorno periódicas
        if (i2>mx-1) { i2=0;dtx=wh;}
        if (i2<0) { i2=mx-1;dtx=-wh;}
        if (j2>my-1) { j2=0;dty=wv;}
        if (j2<0) { j2=my-1;dty=-wv;}

        for (var k=0;k< mat_n_list[i2][j2] ;k++) {

          j=mat_list[i2][j2][k];

          //La partícula no ejerce fuerza sobre sí misma
          if (i==j) continue;

          //Distancia entre cada par de partículas
          dx=(px[j]+dtx)-px[i];
          dy=(py[j]+dty)-py[i];
          d=Math.sqrt(dx*dx+dy*dy);

          ci=pcol[i];cj=pcol[j];
          R1=matR1[ci][cj];
          R2=matR2[ci][cj];
          
          //Ignora las interacciones si la
          //distancia entre partículas es mayor que R2
          if (d==0 || d>R2) continue;
          
          id=K/d;
          
          /*
          //Tipo I (7969ms @ 300)
          ff=0;
          if (d<R1) {
            ff=matcola[ci][cj]*(1-d/R1); 
          } 
          ff+=matcol[ci][cj]*(1-d/R2);
          */
          
          //Tipo 2 (8722ms @ 300)
          if (d<R1) {
            ff=matcola[ci][cj]*(1-d/R1);
          } else {
            Rm=(R1+R2)*0.5;
            R1dR2=R1/R2;
            if (d<Rm) {
              ff=-2*matcol[ci][cj]*(d/R1-1)*(R1dR2)/(1-R1dR2);
            } else {
              ff=2*matcol[ci][cj]*(d/R2-1)/(1-R1dR2);
            }
          }
          
          ffid=ff*id;

          //Acumulación de la interacción
          //sobre la partícula i-ésima
          ax[i]+=dx*ffid;
          ay[i]+=dy*ffid;

        }

      }
    }
    
  }

}

//Generación inicial de las 
//propiedades de las partículas
function init() {

  for (var i=0;i<n;i++) {
    px[i]=Math.random()*wh;
    py[i]=Math.random()*wv;
    vx[i]=0;
    vy[i]=0;
    ax[i]=0;
    ay[i]=0;
    pcol[i]=Math.floor(Math.random()*pt);
  }

  //Genera la matriz de interacción
  //entre tipos de partículas
  matriz_de_interaccion();

}

//Genera la Matriz de interacción
//entre tipos de partículas
function matriz_de_interaccion() {

  //Escala 
  var sc=0.5+Math.random();
  maxR2dR1=-1;

  for (var i=0;i<pt;i++) {

    matcol[i]=[];
    matcola[i]=[];
    matR1[i]=[];
    matR2[i]=[];

    for (var j=0;j<pt;j++) {

      //Constante de interacción i-j
      matcol[i][j]=0.3+0.7*Math.random();
      if (Math.random()<0.5) {
        matcol[i][j]*=-1;
      }

      //Constante de repulsión i-j
      matcola[i][j]=-3*Math.abs(matcol[i][j]);

      //Definición del radio de repulsión (R1) y
      // del radio de interacción (R2)
      matR1[i][j]=(30+40*Math.random())*sc*global_scale;
      matR2[i][j]=(70+180*Math.random())*sc*global_scale;
      if (matR2[i][j]/matR1[i][j]>maxR2dR1) maxR2dR1=matR2[i][j]/matR1[i][j];

    }
  }

}

//Actualiza el tamaño del canvas si se cambia
//el tamaño de la ventana del navegador
function windowResized() {
  wh=windowWidth-10;
  wv=windowHeight-10;
  init_lista();
  resizeCanvas(wh, wv); 
} 

//Definición de las matrices de celdas
//Número y lista de partículas
function init_lista() {
  
  var xi=0;var yi=0;

  for (var i=0;i<mx;i++) {

    mat_list[i]=[];
    mat_n_list[i]=[];
    mat_n_list_zeros[i]=[];

    for (var j=0;j<my;j++) {

      mat_n_list[i][j]=0;
      mat_n_list_zeros[i][j]=0;
      mat_list[i][j]=[];

      for (var k=0;k<400;k++) {
        mat_list[i][j][k]=0;
      }

    }

  }

}
