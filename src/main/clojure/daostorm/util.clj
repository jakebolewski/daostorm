(ns daostorm.util
  (:import [[edu.mcmaster.daostorm ]]))

(def ^:const test-tiff
  "/home/jake/STORM_DATA/comp.tif")

(defn load-imp [filename]
  (.openImage (ij.io.Opener.) filename))

(defn ^ij.gui.PointRoi peak-to-point [peak]
  (let [x (+ (.getXCenter peak) 0.5)
        y (+ (.getYCenter peak) 0.5)]
    (ij.gui.PointRoi. x y)))

(defn ^ij.gui.Overlay peaks-overlay [^java.util.ArrayList peaks]
  (let [overlay (ij.gui.Overlay.)]
    (doseq [p peaks]
      (.add overlay (peak-to-point p)))
    overlay))

(defn ^java.util.ArrayList find-local-maxima
  [img-array taken threshold radius sigma
   background image-width image-height margin]
  (edu.mcmaster.daostorm.Util/findLocalMaxima
    img-array taken threshold radius background
    sigma image-width image-height margin))

(defn imp-to-1d-double-array [^ij.ImagePlus imp]
  (let [img-width (.getWidth imp)
        img-height (.getHeight imp)
        size (* img-width img-height)
        img-array-float (.. imp getProcessor getFloatArray)
        img-array-double (make-array Double/TYPE size)]
    (dotimes [i img-height]
      (dotimes [j img-width]
        (let [idx (+ (* i img-width) j)
              pixel (double (aget img-array-float j i))]
          (aset-double img-array-double idx pixel))))
    img-array-double))

(comment
(defn fconcat [& arrays]
  (let [sizes (map alength arrays)
        sizes-r (vec (reductions + sizes))
        offsets (cons 0 (drop-last sizes-r))
        total (last sizes-r)
        out (make-array (.getComponentType (class (first arrays))) total)]
    (dorun
      (map #(System/arraycopy %1 0 out %2 %3) arrays offsets sizes))
    out))

(defn array? [x]  (-> x class (.isArray)))
(defn see [x] (if (array? x) (map see x) x))

(defn imp-to-1d-double-array2 [^ij.ImagePlus imp]
  (let [img-width (.getWidth imp)
        img-height (.getHeight imp)
        size (* img-width img-height)
        img-array-float (.. imp getProcessor getFloatArray)
        img-array-out (float-array size)]
    (dotimes [i img-height]
      (System/arrayCopy  (aget ^objects img-array-float i) 0
                         img-array-out (* i img-width) img-width))
    img-array-out))
)

(defn show-imp-overlay [^ij.ImagePlus imp
                        ^ij.gui.Overlay overlay]
  (let [_ (ij.ImageJ.)]
    (.setOverlay imp overlay)
    (.show imp)))

(defn local-max-peaks [^ij.ImagePlus imp]
  (let [img-array (imp-to-1d-double-array imp)
        margin 10
        img-width (.getWidth imp)
        img-height (.getHeight imp)
        taken (int-array (alength img-array))
        threshold 300.0
        radius 3.0
        background 100.0
        sigma 1.2]
  (find-local-maxima
        img-array taken threshold radius background
        sigma img-width img-height margin)))

(defn show-local-maxima [filename]
  (let [imp (load-imp filename)
        peaks (local-max-peaks imp)]
    (show-imp-overlay imp (peaks-overlay peaks))))

(comment
  (show-local-maxima test-tiff)
  )

(defn gaussian-blur [imp sigma]
  (let [float-proc (.. imp getProcessor duplicate convertToFloat)
        blured-proc (doto float-proc
                      (.blurGaussian sigma))]
    (doto (ij.ImagePlus.)
      (.setProcessor blured-proc))))

(defn imp-mean [imp]
  (.. imp getStatistics mean))

(defn imp-min [^ij.ImagePlus imp]
  (.. imp getMin))

(defn imp-max [^ij.ImagePlus imp]
  (.. imp getMax))

(defn add-imp [imp1 imp2]
  ;; todo
  (.duplicate imp1))

(defn sub-imp [imp1 imp2]
  ;; todo
  (.duplicate imp2))

(defn add-to-imp [imp x]
  (doto imp (.add x)))

(defn sub-from-imp [imp x]
  (doto imp (.subtract x)))

(defn converged-peaks [peaks]
 ;; todo multi.getConvergedPeaks peaks 0.0 0.0)
  peaks)

(defn get-good-peaks [res threshold sigma]
  res)

(defn remove-close-peaks [peaks sigma neighboorhood]
  peaks)

(defn iterate-fit [fit-func fit-data]
  (let [result (fit-func (:image fit-data))
       fit-peaks (get-good-peaks
                       (get result 0)
                       (* 0.9 (:fit-threshold fit-data))
                       (* 0.5 (:sigma fit-data)))
       ;; remove close peaks
       fit-peaks (remove-close-peaks fit-peaks
                   (:sigma fit-data)
                   (:neighborhood fit-data))
       ;; update fit
       result (fit-func (:image fit-data) fit-peaks)
       fit-peaks (get-good-peaks
                   (get result 0)
                   (:sigma fit-data)
                   (:neighborhood fit-data))]
  (merge fit-data {:peaks fit-peaks
                   :residual (get result 1)})))

(defn update-background-cutoff [fit-data]
  (let [residual (:residual fit-data)
        threshold (:cur-threshold fit-data)
        residual-bg (gaussian-blur residual 8)
        mean-residual-bg (imp-mean residual-bg)
        residual (add-to-imp (sub-imp residual residual-bg)
                              mean-residual-bg)
        background (mean-imp residual)
        cutoff (+ background threshold)]
    (merge fit-data {:residual residual
                     :background background
                     :cutoff cutoff})))

(defn below-threshold [fitdata]
  (let [cur-threshold (:cur-threshold fitdata)
        threshold (:threshold fit-data)]
    (if (> cur-threshold threshold)
      [(- cur-threshold threshold) false]
      [cur-threshold true])))

(defn -main [& args]
  (let [ij (ij.ImageJ.)
        imp (.open (ij.io.Opener.) test-tiff)
        imp2 (.. imp getProcessor convertToFloat)]
    (do (.show imp)
        (.blurGaussian imp2 8)
        (.show (ij.ImagePlus. imp2)))))
