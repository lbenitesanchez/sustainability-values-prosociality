# Draw the estimated M2 structure, without fitting or altering a model.
# Private mode reads an archived fit. Public mode reads the exported parameter
# manifest, model syntax and typed metadata only; no individual data or RDS needed.
# Run: M2_PUBLIC_INPUT_DIR=/path/to/inputs Rscript --vanilla 01_draw_verified_M2.R
# M2_FIGURE_DIR optionally changes the output directory.
suppressPackageStartupMessages(library(lavaan))
library(grid)
model_dir <- Sys.getenv("M2_MODEL_DIR", unset = "/Users/luisbenites/Documents/Codex/Paper con MA - 1/reproduccion-integral-V3-2026-10-10/proyecto/sem-alternativas-2026-10-08/resultados")
public_dir <- Sys.getenv("M2_PUBLIC_INPUT_DIR", unset = "")
out <- Sys.getenv("M2_FIGURE_DIR", unset = "/Users/luisbenites/Documents/Codex/Paper con MA - 1/V3_complete_M2_2026-10-11/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (nzchar(public_dir)) {
  source_mode <- "public exported structure"
  manifest_path <- file.path(public_dir,"M2_parameter_manifest.csv")
  meta_path <- file.path(public_dir,"M2_model_metadata.csv")
  syntax_path <- file.path(public_dir,"M2_puentes2_syntax.txt")
  stopifnot(file.exists(manifest_path),file.exists(meta_path),file.exists(syntax_path))
  pt <- read.csv(manifest_path,stringsAsFactors=FALSE,check.names=FALSE)
  md <- read.csv(meta_path,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(!anyDuplicated(md$field),all(c("sample_n","std.lv","auto.fix.first","fixed.x","meanstructure","estimator","se") %in% md$field))
  converted <- lapply(seq_len(nrow(md)),function(i) {
    if(md$type[i]=="logical") as.logical(md$value[i]) else if(md$type[i]=="integer") as.integer(md$value[i]) else md$value[i]
  })
  names(converted) <- md$field
  opts <- converted
  sample_n <- opts$sample_n
  input_sources <- c(manifest_path,meta_path,syntax_path)
} else {
  source_mode <- "private archived fitted model"
  fit_path <- file.path(model_dir, "M2_puentes2_SEM_MLR.rds")
  syntax_path <- file.path(model_dir, "M2_puentes2_syntax.txt")
  stopifnot(file.exists(fit_path), file.exists(syntax_path))
  fit <- readRDS(fit_path)
  pt <- parameterTable(fit)
  opts <- lavInspect(fit, "options")
  sample_n <- as.integer(lavInspect(fit,"nobs"))
  input_sources <- c(fit_path,syntax_path)
}
stopifnot(sample_n==653L)
# Export the typed, non-record-level metadata needed by public regeneration.
metadata_fields <- c("sample_n","estimator","requested_estimator","se","std.lv","auto.fix.first","fixed.x","meanstructure")
metadata_values <- c(sample_n,opts$estimator,"MLR",opts$se,opts$std.lv,opts$auto.fix.first,opts$fixed.x,opts$meanstructure)
metadata_types <- c("integer","character","character","character",rep("logical",4))
write.csv(data.frame(field=metadata_fields,value=metadata_values,type=metadata_types),file.path(out,"M2_model_metadata.csv"),row.names=FALSE)
if(normalizePath(syntax_path,mustWork=TRUE) != normalizePath(file.path(out,"M2_puentes2_syntax.txt"),mustWork=FALSE)) file.copy(syntax_path,file.path(out,"M2_puentes2_syntax.txt"),overwrite=TRUE)
values <- c("Benevolencia", "Universalismo", "Poder", "Logro")
control <- c("C_Colombia", "C_Peru", "C_Mexico", "Sexo_num", "Estudios_num", "Edad_z", "NSE_z")
endogenous <- c(values, "SC", "SBL", "PSB")

# Verify every drawn or captioned relation against the saved fit, not a sketch.
expected_focal <- rbind(
  expand.grid(lhs=values, rhs="CSR1_bin"),
  expand.grid(lhs=c("SC","SBL"), rhs=c(values,"CSR1_bin")),
  expand.grid(lhs="PSB", rhs=c("SC","SBL",values,"CSR1_bin")))
actual_focal <- pt[pt$op == "~" & !pt$rhs %in% control, c("lhs","rhs")]
key <- function(z) paste(z$lhs, z$rhs, sep="|")
stopifnot(setequal(key(expected_focal),key(actual_focal)), nrow(actual_focal)==21L)
syntax_pt <- lavaanify(paste(readLines(syntax_path, warn=FALSE),collapse="\n"),model.type="sem",fixed.x=TRUE)
stopifnot(setequal(key(syntax_pt[syntax_pt$op=="~" & !syntax_pt$rhs %in% control,]),key(actual_focal)))
actual_control <- pt[pt$op == "~" & pt$rhs %in% control,c("lhs","rhs")]
expected_control <- expand.grid(lhs=endogenous,rhs=control)
stopifnot(setequal(key(expected_control),key(actual_control)),nrow(actual_control)==49L)
expected_measurement <- rbind(
  data.frame(lhs="SC",rhs=c("APS1","APS2","APS3","APS4","APS5","APS6","APS7","APS9","APS10","APS14","APS13")),
  data.frame(lhs="SBL",rhs=c("APS11","APS12","APS13","APS15","APS16","APS18","APS19","APS20","APS21","APS2")),
  data.frame(lhs="PSB",rhs=paste0("PP",1:3)))
stopifnot(setequal(key(expected_measurement),key(pt[pt$op=="=~",])),sum(pt$op=="=~")==24L)
stopifnot(setequal(key(syntax_pt[syntax_pt$op=="=~",]),key(pt[pt$op=="=~",])))
expected_cov <- c(apply(combn(values,2),2,paste,collapse="|"),"SC|SBL")
actual_cov <- pt[pt$op=="~~" & pt$lhs!=pt$rhs & pt$free>0,c("lhs","rhs")]
stopifnot(setequal(expected_cov,key(actual_cov)),nrow(actual_cov)==7L)
markers <- pt[pt$op=="=~" & pt$free==0,c("lhs","rhs","ustart")]
stopifnot(setequal(key(markers),c("SC|APS1","SBL|APS11","PSB|PP1")),all(markers$ustart==1),!opts$std.lv,opts$fixed.x)

# Manifest includes *all* parameters, including fixed exogenous sample moments,
# intercepts, residual variances and defined indirect products not drawn as edges.
manifest <- pt
manifest$representation <- "not drawn; see model specification"
manifest$representation[pt$op=="~" & !pt$rhs %in% control] <- "panel A: regression arrow"
manifest$representation[pt$op=="~" & pt$rhs %in% control] <- "caption: seven controls in every endogenous equation"
manifest$representation[pt$op=="=~"] <- "panel B: loading arrow, grouped where stated"
manifest$representation[pt$op=="~~" & pt$lhs!=pt$rhs & pt$free>0] <- "caption: residual covariance"
manifest$representation[pt$op=="~~" & pt$lhs==pt$rhs] <- "caption: variances/residual variances omitted"
manifest$representation[pt$op=="~1"] <- "caption: intercepts omitted"
manifest$representation[pt$op==":="] <- "derived product; not an additional estimated arrow"
manifest$representation[pt$op=="~~" & pt$free==0] <- "fixed exogenous moment; not a free covariance arrow"
write.csv(manifest,file.path(out,"M2_parameter_manifest.csv"),row.names=FALSE,na="")
checks <- data.frame(check=c("sample_n","focal_regressions","control_regressions","measurement_loadings","free_residual_covariances","marker_loadings_fixed_to_1","std_lv","fixed_x","model_refit"),
  value=c(sample_n,nrow(actual_focal),nrow(actual_control),sum(pt$op=="=~"),nrow(actual_cov),nrow(markers),opts$std.lv,opts$fixed.x,"NO"))
write.csv(checks,file.path(out,"M2_figure_checks.csv"),row.names=FALSE)
write.csv(data.frame(file=basename(input_sources),md5=unname(tools::md5sum(input_sources))),file.path(out,"M2_input_hashes.csv"),row.names=FALSE)

# Portable SVG writer records the same R/grid coordinates used by the PDF.
# This avoids a Cairo/XQuartz dependency on macOS and preserves vector geometry.
svg_capture <- FALSE
svg_parts <- character()
xml_escape <- function(x) { x <- gsub("&", "&amp;", x, fixed=TRUE); x <- gsub("<", "&lt;", x, fixed=TRUE); gsub(">", "&gt;", x, fixed=TRUE) }
xy_svg <- function(x,y) {
  p <- deviceLoc(unit(x,"npc"),unit(y,"npc"),valueOnly=TRUE)
  data.frame(x=p$x*72,y=(9.5-p$y)*72)
}
add_svg <- function(x) if (svg_capture) svg_parts <<- c(svg_parts,x)
text_draw <- function(label,x,y,just="centre",gp=gpar()) {
  grid.text(label,x,y,just=just,gp=gp)
  if(svg_capture) {
    p<-xy_svg(x,y); fs<-if(is.null(gp$fontsize))12 else gp$fontsize
    color<-if(is.null(gp$col))"black" else gp$col
    weight<-if(!is.null(gp$font)&&gp$font==2)"bold" else "normal"
    anchor<-if(identical(just,"left"))"start" else if(identical(just,"right"))"end" else "middle"
    ls<-strsplit(as.character(label),"\n",fixed=TRUE)[[1]]
    for(i in seq_along(ls)) add_svg(sprintf('<text x="%.3f" y="%.3f" text-anchor="%s" dominant-baseline="central" font-family="Arial,Helvetica,sans-serif" font-size="%.3f" font-weight="%s" fill="%s">%s</text>',p$x,p$y+(i-(length(ls)+1)/2)*fs*1.05,anchor,fs,weight,color,xml_escape(ls[i])))
  }
}
rect_draw <- function(x,y,width,height,r,gp) {
  grid.roundrect(x,y,width,height,r=r,gp=gp)
  if(svg_capture) {
    p<-xy_svg(x-width/2,y+height/2); q<-xy_svg(x+width/2,y-height/2)
    add_svg(sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" rx="2" fill="%s" stroke="%s" stroke-width="%.2f"/>',p$x,p$y,q$x-p$x,q$y-p$y,gp$fill,gp$col,gp$lwd*.75))
  }
}
polygon_draw <- function(x,y,gp) {
  grid.polygon(x,y,gp=gp)
  if(svg_capture) {
    p<-xy_svg(x,y)
    add_svg(sprintf('<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.2f"/>',paste(sprintf('%.3f,%.3f',p$x,p$y),collapse=' '),gp$fill,gp$col,gp$lwd*.75))
  }
}

# Drawing primitives use normalized coordinates within each panel viewport.
ink <- "#172B3A"; blue <- "#1C6588"; orange <- "#AD5A2D"; gray <- "#65727B"; teal <- "#1D736F"
line <- function(x,y,col=ink,lwd=1.5,lty=1,ends="last") {
  grid.lines(x,y,default.units="npc",gp=gpar(col=col,lwd=lwd,lty=lty),arrow=arrow(length=unit(.085,"inches"),ends=ends,type="closed"))
  if(svg_capture) {
    p<-xy_svg(x,y); dash<-if(lty==2)' stroke-dasharray="5 4"' else if(lty==3)' stroke-dasharray="1.5 3"' else ''
    add_svg(sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%.2f"%s marker-end="url(#a%s)"/>',paste(sprintf('%.3f,%.3f',p$x,p$y),collapse=' '),col,lwd*.75,dash,substring(col,2)))
  }
}
curve <- function(x,y,col=ink,lwd=1.5,lty=1){
  t<-seq(0,1,length.out=101); b<-cbind((1-t)^3,3*(1-t)^2*t,3*(1-t)*t^2,t^3)
  line(drop(b%*%x),drop(b%*%y),col,lwd,lty)
}
rect_node <- function(x,y,label,w=.18,h=.085,fill="white",size=13.5){
  rect_draw(x,y,w,h,r=unit(.025,"inches"),gp=gpar(fill=fill,col=ink,lwd=1.35))
  text_draw(label,x,y,gp=gpar(fontsize=size,col=ink,fontfamily="sans"))
}
latent <- function(x,y,label,w=.22,h=.11,size=13.5){
  theta<-seq(0,2*pi,length.out=100)
  polygon_draw(x+w/2*cos(theta),y+h/2*sin(theta),gp=gpar(fill="#EDF5F6",col=ink,lwd=1.5))
  text_draw(label,x,y,gp=gpar(fontsize=size,col=ink,fontfamily="sans"))
}

draw <- function(){
  grid.newpage()
  text_draw("M2: estimated structural and measurement specification",x=.035,y=.973,just="left",gp=gpar(fontsize=18,fontface="bold",col=ink))
  text_draw("Rectangles: observed variables    Ovals: latent factors    No coefficient estimates are shown",x=.035,y=.946,just="left",gp=gpar(fontsize=11.5,col=gray))
  pushViewport(viewport(x=.5,y=.665,width=.93,height=.50))
  text_draw("A   Structural regressions",x=0,y=1.02,just="left",gp=gpar(fontsize=14,fontface="bold",col=ink))
  xs<-c(.18,.39,.61,.82)
  # CSR1 -> each observed value.
  for(x in xs) line(c(.5,x),c(.885,.715),gray,1.5)
  # Each value -> each attitude factor: eight individual, unpooled paths.
  for(i in seq_along(xs)) {
    off <- c(-.05,-.015,.015,.05)[i]
    line(c(xs[i],.36 + off),c(.625,.365+.065*sqrt(1-(off/.115)^2)),blue,1.5)
    line(c(xs[i],.64 + off),c(.625,.365+.080*sqrt(1-(off/.1225)^2)),blue,1.5)
  }
  # Direct value -> PSB paths, routed around attitude nodes for legibility.
  curve(c(xs[1]-.055,.025,.02,.425),c(.625,.26,.05,.131),orange,1.4,2)
  curve(c(xs[2]-.04,.11,.14,.465),c(.625,.30,.16,.142),orange,1.4,2)
  curve(c(xs[3]+.04,.89,.86,.535),c(.625,.30,.16,.142),orange,1.4,2)
  curve(c(xs[4]+.055,.975,.98,.575),c(.625,.26,.05,.131),orange,1.4,2)
  # CSR1 -> SC, SBL and PSB (three distinct fitted direct regressions).
  curve(c(.4,-.03,-.01,.245),c(.925,.98,.39,.365),gray,1.4,3)
  curve(c(.6,1.03,1.01,.7625),c(.925,.98,.39,.365),gray,1.4,3)
  curve(c(.59,1.10,1.10,.622),c(.955,1.08,-.03,.075),gray,1.4,3)
  # SC and SBL -> PSB.
  line(c(.36,.47),c(.30,.145),teal,2)
  line(c(.64,.53),c(.285,.145),teal,2)
  rect_node(.5,.93,"CSR1\nCourse attendance",.20,.09,fill="#F4F5F6",size=13)
  labs<-c("Benevolence","Universalism","Power","Achievement")
  for(i in seq_along(xs)) rect_node(xs[i],.67,labs[i],.18,.09,fill="#FCF8F0")
  latent(.36,.365,"SC\nSystemic consciousness",.23,.13,size=12.2)
  latent(.64,.365,"SBL\nSustainable business\nleadership",.245,.16,size=11.8)
  latent(.5,.08,"PSB\nProsocial self-descriptions",.245,.13,size=12.1)
  popViewport()

  # Compact style key. Colors distinguish arrow families, not significance.
  pushViewport(viewport(x=.5,y=.383,width=.93,height=.025))
  line(c(.01,.055),c(.5,.5),blue);text_draw("Values -> SC/SBL",.065,.5,just="left",gp=gpar(fontsize=10.6,col=ink))
  line(c(.30,.345),c(.5,.5),orange,lty=2);text_draw("Values -> PSB",.355,.5,just="left",gp=gpar(fontsize=10.6,col=ink))
  line(c(.555,.60),c(.5,.5),teal,lwd=2);text_draw("SC/SBL -> PSB",.61,.5,just="left",gp=gpar(fontsize=10.6,col=ink))
  line(c(.82,.865),c(.5,.5),gray,lty=3);text_draw("CSR1 direct paths",.875,.5,just="left",gp=gpar(fontsize=10.2,col=ink))
  popViewport()

  pushViewport(viewport(x=.5,y=.218,width=.93,height=.265))
  text_draw("B   Measurement (the same SC, SBL and PSB factors as in panel A)",0,1.02,just="left",gp=gpar(fontsize=14,fontface="bold",col=ink))
  # Primary loading arrows. Grouped rectangles enumerate individual indicators.
  line(c(.15,.29),c(.73,.845),blue)
  line(c(.15,.33),c(.73,.595),blue)
  line(c(.15,.29),c(.27,.155),blue)
  line(c(.15,.33),c(.27,.405),blue)
  # Exact two added cross-loadings of M2.
  line(c(.15,.33),c(.715,.42),orange,lty=2)
  line(c(.15,.33),c(.285,.58),orange,lty=2)
  latent(.10,.75,"SC",.10,.18,size=13.5)
  latent(.10,.25,"SBL",.10,.18,size=13.5)
  rect_node(.45,.845,"APS1, APS3, APS4, APS5, APS6,\nAPS7, APS9, APS10, APS14",.32,.20,size=11.1)
  rect_node(.40,.595,"APS2",.14,.13,size=12)
  rect_node(.40,.405,"APS13",.14,.13,size=12)
  rect_node(.45,.155,"APS11, APS12, APS15, APS16,\nAPS18, APS19, APS20, APS21",.32,.20,size=11.1)
  text_draw("Primary loadings",.05,.51,just="left",gp=gpar(fontsize=10,col=blue))
  text_draw("Dashed orange: cross-loadings\nSC -> APS13; SBL -> APS2",.36,-.015,just="centre",gp=gpar(fontsize=10.6,col=orange))
  # Prosocial measurement: three observed responses, not three SEM outcomes.
  latent(.74,.50,"PSB",.105,.20,size=13.5)
  for(y in c(.8,.5,.2)) line(c(.79,.88),c(.5,y),blue)
  rect_node(.935,.8,"PP1\nHelping",.12,.18,size=11.1)
  rect_node(.935,.5,"PP2\nSharing",.12,.18,size=11.1)
  rect_node(.935,.2,"PP3\nEmpathic support",.12,.18,size=10.5)
  popViewport()
  text_draw("Grouped boxes abbreviate separate observed indicators; they are not composite scores. All relations are cross-sectional.",.035,.035,just="left",gp=gpar(fontsize=11.2,col=gray))
}

# Quartz creates the correct pixels but can write a 72 dpi metadata tag.
# Write the pHYs PNG chunk explicitly; pixels are unchanged. CRC is checked by PNG readers.
set_png_dpi <- function(path,dpi) {
  z<-readBin(path,"raw",n=file.info(path)$size)
  integer_raw<-function(n) writeBin(as.integer(n),raw(),size=4,endian="big")
  ppm<-round(dpi/.0254)
  dat<-c(integer_raw(ppm),integer_raw(ppm),as.raw(1))
  chunk_body<-c(charToRaw("pHYs"),dat)
  crc_hex<-digest::digest(chunk_body,algo="crc32",serialize=FALSE)
  crc<-as.raw(strtoi(substring(crc_hex,seq(1,7,2),seq(2,8,2)),16L))
  new_chunk<-c(integer_raw(9),chunk_body,crc)
  keep<-z[1:8]; pos<-9L
  while(pos<=length(z)) {
    len<-readBin(z[pos:(pos+3L)],integer(),n=1,size=4,endian="big")
    end<-pos+len+11L
    typ<-rawToChar(z[(pos+4L):(pos+7L)])
    if(typ!="pHYs") keep<-c(keep,z[pos:end])
    if(typ=="IHDR") keep<-c(keep,new_chunk)
    pos<-end+1L
  }
  writeBin(keep,path)
}

# Both vector formats use the same coordinates. PNG uses the native macOS device
# when available; Linux/Windows use Cairo. No external editor modifies the image.
svg_capture <- TRUE
pdf(file.path(out,"Figure_M2_verified.pdf"),width=13,height=9.5,useDingbats=FALSE,family="Helvetica");draw();dev.off()
svg_capture <- FALSE
svg_defs <- vapply(c(ink,blue,orange,gray,teal),function(co)sprintf('<marker id="a%s" markerWidth="6" markerHeight="6" refX="5.4" refY="3" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0 L6,3 L0,6 Z" fill="%s"/></marker>',substring(co,2),co),character(1))
writeLines(c('<svg xmlns="http://www.w3.org/2000/svg" width="13in" height="9.5in" viewBox="0 0 936 684">','<rect width="936" height="684" fill="white"/>','<defs>',svg_defs,'</defs>',svg_parts,'</svg>'),file.path(out,"Figure_M2_verified.svg"))
bitmap_type <- if(Sys.info()[["sysname"]]=="Darwin") "quartz" else "cairo"
png(file.path(out,"Figure_M2_verified_600dpi.png"),width=13,height=9.5,units="in",res=600,type=bitmap_type);draw();dev.off()
set_png_dpi(file.path(out,"Figure_M2_verified_600dpi.png"),600)
png(file.path(out,"Figure_M2_preview.png"),width=1950,height=1425,res=150,type=bitmap_type);draw();dev.off()
set_png_dpi(file.path(out,"Figure_M2_preview.png"),150)

caption <- paste(
"Figure 1. M2 structural regressions and measurement specification.",
"Panel A shows all 21 focal regressions; arrow styles distinguish predictor families, not statistical significance.",
"The four value portraits are separate observed variables, not indicators of an additional latent factor.",
"Panel B shows the measurement specification: SC and SBL have 19 distinct APS-S indicators and the two cross-loadings SC → APS13 and SBL → APS2; PSB is measured by PP1–PP3.",
"Each grouped box abbreviates separate indicators and loadings, not a score or parcel.",
"Seven covariates (Colombia, Peru and Mexico institutional indicators, sex, study level, standardized age and standardized socioeconomic level) enter every equation for the four values, SC, SBL and PSB (49 control regressions, omitted for readability).",
"The six residual covariances among the four values and the SC–SBL disturbance covariance are freely estimated; these are conditional residual covariances, not total factor correlations.",
"Observed residual variances, endogenous disturbances, intercepts and fixed exogenous sample moments are omitted. APS1, APS11 and PP1 are the marker loadings fixed to 1 for SC, SBL and PSB, respectively.",
"The diagram follows the archived M2-SEM estimated with MLR (N = 653), not the model treating PP1–PP3 as three outcomes; complete loading and path estimates are reported in the supplement.",
"Arrows denote cross-sectional regressions and do not establish causality or temporal order.")
writeLines(caption,file.path(out,"Figure_M2_caption.txt"))
caption_short <- paste(
"Figure 1. M2 structural and measurement specification (MLR, N = 653).",
"Panel A shows the 21 focal regressions; panel B shows measurement, with SC → APS13 and SBL → APS2 as the two cross-loadings.",
"Rectangles denote observed variables and ovals denote latent factors; the four value portraits are separate observed variables. Grouped boxes abbreviate individual indicators, not composites.",
"Seven controls enter each equation for the four values, SC, SBL and PSB: three institutional indicators, sex, study level, standardized age and socioeconomic level.",
"Six residual covariances among values and the SC–SBL disturbance covariance are estimated but omitted; residual variances, intercepts and fixed exogenous moments are also omitted.",
"APS1, APS11 and PP1 are marker loadings fixed to 1. Colors identify arrow families, not significance. Arrows denote cross-sectional associations, not causality or temporal order.")
writeLines(caption_short,file.path(out,"Figure_M2_caption_short.txt"))
writeLines(c(
  "Resumen de la figura M2", "",
  paste("Origen de la estructura:",source_mode,". Se contrastó con la sintaxis archivada; no se reestimó ningún modelo. El modo público regenera la figura desde exportaciones y no reproduce la estimación desde datos individuales."),
  "El panel A representa las 21 regresiones focales; el panel B representa 24 cargas: 22 principales más las dos cruzadas de M2.",
  "Los cuatro retratos de valores aparecen como rectángulos observados separados. SC, SBL y PSB aparecen como óvalos. PP1–PP3 son indicadores observados de PSB en este SEM.",
  "La leyenda documenta las 49 regresiones de control, siete covarianzas residuales libres, las tres cargas marcadoras fijadas en 1 y los elementos omitidos por legibilidad.",
  "La clasificación completa de todos los parámetros está en M2_parameter_manifest.csv; los recuentos comprobados están en M2_figure_checks.csv.",
  "Los archivos PDF y SVG son vectoriales; el PNG de entrega mide 7800 × 5700 píxeles a 600 ppp. Figure_M2_preview.png es sólo una vista previa.",
  "No se utilizan índices, coeficientes ni resultados de otro modelo. Los colores y estilos no codifican significación estadística."
),file.path(out,"RESUMEN_FIGURA_M2_ES.txt"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo_figure.txt"))
cat("Verified M2 figure complete. Input mode:",source_mode,". No model fitted.\n")
