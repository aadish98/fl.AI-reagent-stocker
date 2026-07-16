const pptxgen = require("pptxgenjs");
let pres = new pptxgen();
let slide = pres.addSlide();
slide.addText("Hello World!", { x: 0.5, y: 0.5, fontSize: 36, color: "363636" });
pres.writeFile({ fileName: "TestMin.pptx" }).then(() => console.log("done"));
