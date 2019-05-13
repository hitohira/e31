var gl;
var programInfo;

//extern var fits;
//extern var color_f;
//extern var pres;
//extern var color_p;

var fit_en = true;
var pre_en = true;
var down = false;
var down_x;
var down_y;
/////////////////////////////////////////
/// mouse event
/////////////////////////////////////////
function checkCurve(){
	fit_en = document.form.curve[0].checked;
	pre_en = document.form.curve[1].checked;
	drawScene();
}

function mouseDown(event){
	down = true;
	down_x = getX(event);
	down_y = getY(event);
}

function mouseUp(){
	var x = getX(event);
	var y = getY(event);
	if(y - down_y > 0.5){ // 拡大
		expand(1.5);
	}
	else if(y - down_y > 0.2){ // 拡大
		expand(1.2);
	}
	else if(y - down_y < -0.5){ //縮小
		expand(1.0/1.5);
	}
	else if(y - down_y < -0.2){ // 拡大
		expand(1.0/1.2);
	}
	console.log("hello");
	down = false;

	drawScene();
}

function dblClick(event){ // centerをここにする。
	var x = getX(event);
	var y = getY(event);
	
	for(var i = 0; i < fits.length/3; i++){
		fits[i*3] -= x;
		fits[i*3+1] -= y;
	}
	for(var i = 0; i < pres.length/3; i++){
		pres[i*3] -= x;
		pres[i*3+1] -= y;
	}
	drawScene();
}

function getX(event){
		const x = event.clientX;
			return (1.0 * x / gl.canvas.clientWidth - 0.5) * 2;
}
function getY(event){
		const y = event.clientY; 	
			return (1.0 * y / gl.canvas.clientHeight - 0.5) * -2
}

function expand(times){
	for(var i = 0; i < fits.length/3; i++){
		fits[i*3] *= times;
		fits[i*3+1] *= times;
	}
	for(var i = 0; i < pres.length/3; i++){
		pres[i*3] *= times;
		pres[i*3+1] *= times;
	}
	
}
/////////////////////////////////////////
/// main 
/////////////////////////////////////////

function main(){
	initCanvas(); // init gl
	const vsSource = vsCode();
	const fsSource = fsCode();
	const shaderProgram = initShaderProgram(vsSource,fsSource); // init shaderProgram
	programInfo = {
		program: shaderProgram,
		attribLocations: {
			vertexPosition: gl.getAttribLocation(shaderProgram,"aVertexPosition"),
			vertexColor: gl.getAttribLocation(shaderProgram,"aVertexColor"),
		},
	};

	drawScene();

}

function vsCode(){
	return `
		attribute vec4 aVertexPosition;
		attribute vec4 aVertexColor;

		varying lowp vec4 vColor;
		
		void main(){
			gl_Position = aVertexPosition;
			gl_PointSize = 3.0;
			vColor = aVertexColor;
		}
	`;
}

function fsCode(){
	return `
		varying lowp vec4 vColor;

		void main(){
			gl_FragColor = vColor;
		}
	`;
}

function drawScene(){
	gl.clearColor(0.0,0.0,0.0,1.0);
	gl.clearDepth(1.0);
	gl.enable(gl.DEPTH_TEST);
	gl.depthFunc(gl.LEQUAL);

	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	
	
	var fitData = {
		buffers : initBuffers(fits,color_f),
		vCount : fits.length / 3,
	};
	var preData = {
		buffers : initBuffers(pres,color_p),
		vCount : pres.length / 3,
	};
	// draw curve
	if(pre_en){
		{	
			const numComponents = 3;
			const type = gl.FLOAT;
			const normalize = false;
			const stride = 0;
			const offset = 0;
			gl.bindBuffer(gl.ARRAY_BUFFER,preData.buffers.position);
			gl.vertexAttribPointer(programInfo.attribLocations.vertexPosition,
		                       numComponents,
													 type,
													 normalize,
													 stride,
													 offset);
			gl.enableVertexAttribArray(programInfo.attribLocations.vertexPosition);	
		}
		{	
			const numComponents = 4;
			const type = gl.FLOAT;
			const normalize = false;
			const stride = 0;
			const offset = 0;
			gl.bindBuffer(gl.ARRAY_BUFFER,preData.buffers.color);
			gl.vertexAttribPointer(programInfo.attribLocations.vertexColor,
		                       numComponents,
													 type,
													 normalize,
													 stride,
													 offset);
			gl.enableVertexAttribArray(programInfo.attribLocations.vertexColor);	
		}
		gl.useProgram(programInfo.program);
		{
			const offset = 0;
			const vertexCount = preData.vCount;
			gl.drawArrays(gl.LINE_STRIP,offset,vertexCount);
			gl.drawArrays(gl.POINTS,offset,vertexCount);
		}
	}
	if(fit_en){
		{	
			const numComponents = 3;
			const type = gl.FLOAT;
			const normalize = false;
			const stride = 0;
			const offset = 0;
			gl.bindBuffer(gl.ARRAY_BUFFER,fitData.buffers.position);
			gl.vertexAttribPointer(programInfo.attribLocations.vertexPosition,
		                       numComponents,
													 type,
													 normalize,
													 stride,
													 offset);
			gl.enableVertexAttribArray(programInfo.attribLocations.vertexPosition);	
		}
		{	
			const numComponents = 4;
			const type = gl.FLOAT;
			const normalize = false;
			const stride = 0;
			const offset = 0;
			gl.bindBuffer(gl.ARRAY_BUFFER,fitData.buffers.color);
			gl.vertexAttribPointer(programInfo.attribLocations.vertexColor,
		                       numComponents,
													 type,
													 normalize,
													 stride,
													 offset);
			gl.enableVertexAttribArray(programInfo.attribLocations.vertexColor);	
		}
		gl.useProgram(programInfo.program);
		{
			const offset = 0;
			const vertexCount = fitData.vCount;
			gl.drawArrays(gl.LINE_STRIP,offset,vertexCount);
	//		gl.drawArrays(gl.POINTS,offset,vertexCount);
		}
	}
}



/////////////////////////////////////////
/// util functions 
/////////////////////////////////////////

function initBuffers(positions,colors){
	const positionBuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER,positionBuffer);
	gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(positions),gl.STATIC_DRAW);

	const colorBuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);
	gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(colors),gl.STATIC_DRAW);

	return {
		position : positionBuffer,
		color : colorBuffer,
	};
}

function initCanvas(){
	const canvas = document.querySelector("#glcanvas");
	gl = canvas.getContext("webgl");
	if(!gl){
		alert("may not support WebGL");
	}
}


function initShaderProgram(vsSource,fsSource){
	var vs = loadShader(gl.VERTEX_SHADER,vsSource);
	var fs = loadShader(gl.FRAGMENT_SHADER,fsSource);

	const shaderProgram = gl.createProgram();
	gl.attachShader(shaderProgram,vs);
	gl.attachShader(shaderProgram,fs);
	gl.linkProgram(shaderProgram);

	if(!gl.getProgramParameter(shaderProgram,gl.LINK_STATUS)){
		alert("fail at initiarizing shader program");
	}

	return shaderProgram;
}
// idが`id`である<script>からシェーダプログラムを取ってくる
function loadShader(type,source){
	const shader = gl.createShader(type);
	gl.shaderSource(shader,source);
	
	gl.compileShader(shader);
	if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS)){
		alert("fail at compiling shader: " + gl.getShaderInfoLog(shader));
		gl.deleteShader(shader);
		return null;
	}
	return shader;
}

