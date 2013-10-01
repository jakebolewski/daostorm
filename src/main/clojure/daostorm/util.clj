(ns daostorm.util
  (:import [edu.mcmaster.daostorm Util MultiFit PeakStatus]))

(def ^:const test-tiff
  "/home/jake/STORM_DATA/comp.tif")

(def colors {:green java.awt.Color/GREEN
             :red java.awt.Color/RED
             :yellow java.awt.Color/YELLOW
             :blue java.awt.Color/BLACK
             :gray java.awt.Color/LIGHT_GRAY})

(defn load-imp [filename]
  (.openImage (ij.io.Opener.) filename))

(defn ^ij.gui.PointRoi peak->point
  [peak & {:keys [color] :or {:color (:light-gray colors)}}]
  (let [x (+ (.getXCenter peak) 0.5)
        y (+ (.getYCenter peak) 0.5)]
    (doto (ij.gui.PointRoi. x y)
      (.setStrokeColor color))))

(defn ^ij.gui.Overlay peaks-overlay [^java.util.ArrayList peaks]
  (let [overlay (ij.gui.Overlay.)]
    (doseq [p peaks]
      (.add overlay (peak->point p)))
    overlay))

(defn ^java.util.ArrayList find-local-maxima
  [img-array taken threshold radius background sigma
   image-width image-height margin]

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
        threshold 200.0
        radius 1.5
        background 50.0
        sigma 1.5]
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
  (.. imp getStatistics min))

(defn imp-max [^ij.ImagePlus imp]
  (.. imp getStatistics max))

(defn add-imp [imp1 imp2]
  (Util/addImageProcessors
    (.getProcessor imp1) (.getProcessor imp2)))

(defn sub-imp [imp1 imp2]
  (Util/subtractImageProcessors
    (.getProcessor imp1) (.getProcessor imp2)))

(defn add-to-imp [imp x]
  (doto imp (.add x)))

(defn sub-from-imp [imp x]
  (.. imp getProcessor (subtract (double x)))
  imp)

(defn converged-peaks [multi peaks]
  (.getConvergedPeaks multi peaks 0.0 0.0))

(defn get-good-peaks [peaks threshold sigma]
  (Util/getGoodPeaks peaks threshold sigma))

(defn remove-close-peaks [peaks sigma neighboorhood]
  (Util/removeClosePeaks peaks sigma neighboorhood))

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
        background (imp-mean residual)
        cutoff (+ background threshold)]
    (merge fit-data {:residual residual
                     :background background
                     :cutoff cutoff})))

(defn below-threshold [fit-data]
  (let [cur-threshold (:cur-threshold fit-data)
        threshold (:threshold fit-data)]
    (if (> cur-threshold threshold)
      [(- cur-threshold threshold) false]
      [cur-threshold true])))

(defn find-peaks [fit-data]
  ;;todo: support masking with ImageJ's ROI's
  (let [masked-residual (:residual fit-data)]

    ))

(defn do-fit [imp peaks tolerance max-iters z-fit?]
  (let [num-peaks (.size peaks)
        data (imp-to-1d-double-array imp)
        width (.getWidth imp)
        height (.getHeight imp)
        multi (MultiFit. data peaks tolerance width height z-fit?)
        num-iters (atom 0)]
    (do
      (.iter3D multi)
      (loop [i 1]
        (when (and (> (.getNumUnconverged multi))
                   (< i max-iters))
          (.iter3D multi)
          (swap! num-iters inc)
          ;;(prn (str "Total Error: iteration " i " error: " (.getTotalError multi)))
          (recur (inc i))))
      (if (== @num-iters (dec max-iters))
        (prn (str "Failed to converge in: " @num-iters " " (.getNumUnconverged multi)))
        (prn (str "Multi-fit converged in: " @num-iters " " (.getNumUnconverged multi))))
      [(.getResults multi)
       (.getForeground multi)
       (.getResidual multi)])))


(defn peak->oval
  [peak & {:keys [color] :or {color (:gray colors)}}]
  (let [height (* (.getYWidth peak) 2.5)
        width (* (.getXWidth peak) 2.5)
        x (- (+ (.getXCenter peak) 0.5) (/ width 2.0))
        y (- (+ (.getYCenter peak) 0.5) (/ height 2.0))]
    (doto (ij.gui.OvalRoi. x y width height)
      (.setStrokeColor color))))

(defn filter-converged [peaks]
  (filter #(.hasStatus % PeakStatus/CONVERGED) peaks))

(defn converged? [peak]
  (.hasStatus peak PeakStatus/CONVERGED))

(defn ^ij.gui.Overlay peaks-update-overlay [^java.util.ArrayList maxima
                                            ^java.util.ArrayList peaks]
  (let [overlay (ij.gui.Overlay.)]
    (doseq [p maxima]
      (.add overlay (peak->point p :color (:gray colors))))
    (doseq [p peaks]
      (if (converged? p)
        (.add overlay (peak->oval p :color (:green colors)))
        (.add overlay (peak->oval p :color (:red colors)))))
    overlay))

(defn visualize-double-buffer [^doubles residual width height title]
  (let [float-buffer (float-array (alength residual))
        float-proc (ij.process.FloatProcessor. width height)
        imp (doto (ij.ImagePlus.)
              (.setTitle title))]
    (dotimes [i (alength residual)]
      (aset-float float-buffer i (aget residual i)))
    (.setPixels float-proc float-buffer)
    (.setProcessor imp float-proc)
    (.show imp)))

(defn show-local-maxima [filename]
  (let [imp (load-imp filename)
        peaks (local-max-peaks imp)]
    (show-imp-overlay imp (peaks-overlay peaks))))

(defn fit-stats [peaks]
  (let [stats (Util/fitStats peaks)]
    {:converged (aget stats 0)
     :bad (aget stats 1)
     :unconverged (aget stats 2)
     :total (aget stats 3)}))

(defn subtract-baseline [imp baseline]
  {:pre [(>= baseline 0)]}
  (sub-from-imp imp baseline))

(defn subtract-baseline! [imp baseline]
  {:pre [(> baseline 0)]}
  (let [proc (.. imp getProcessor)
        npixels (.getPixelCount proc)]
  (loop [idx npixels]
    (let [val (-> (.getf proc idx)
                  (- baseline)
                  (max 0.0))]
      (.setf proc idx val)))))

(defn test-fitting []
  (let [imp (-> (doto (load-imp test-tiff) (.setSlice 2))
                (subtract-baseline 100))
        width (.getWidth imp)
        height (.getHeight imp)
        peaks (local-max-peaks imp)
        local-max (Util/copyPeakList peaks)
        [result-peaks model residual] (do-fit imp peaks 1e-6 500 false)
        fit-stats (fit-stats result-peaks)]
    (doseq [p (take 100 result-peaks)]
      (prn [(.getXWidth p) (.getYWidth p)]))
    (prn fit-stats)
    (show-imp-overlay imp (peaks-update-overlay local-max result-peaks))
    (visualize-double-buffer residual width height "Residual")
    (visualize-double-buffer model width height "Model")
  ))

;; TODO: raise better error than null pointer when fitting z
;; NullPointerException   edu.mcmaster.daostorm.MultiFit.calcWidthsFromZ (MultiFit.java:227)

(defn test-multi []
  (let [imp (-> (doto (load-imp test-tiff) (.setSlice 2))
                (subtract-baseline 100))
        peaks (local-max-peaks imp)
        local-max (Util/copyPeakList peaks)
        data (imp-to-1d-double-array imp)
        width (.getWidth imp)
        height (.getHeight imp)
        multi (MultiFit. data peaks 1e-6 width height false)
        result (.getResults multi)
        foreground (.getForeground multi)
        residual (.getResidual multi)]
    (show-imp-overlay imp (peaks-update-overlay local-max result))
    (visualize-double-buffer residual width height "Residual")
    (visualize-double-buffer foreground width height "Model")
    ))

(defn gaussian [height xc yc sigma]
  (fn [x y] (* height (Math/exp (- (/ (+ (Math/pow (- x xc) 2) (Math/pow (- y yc) 2))
                                   (* 2.0 (Math/pow sigma 2))))))))

(defn centered-spot [size height sigma]
  (let [proc (ij.process.FloatProcessor. size size)
        gauss-f (gaussian height (/ size 2.0) (/ size 2.0) sigma)]
    (loop [y 0]
      (when (< y size)
        (loop [x 0]
          (when (< x size)
            (.setf proc x y (gauss-f x y))
            (recur (inc x))))
        (recur (inc y))))
    proc))

(defn test-single-centered []
  (let [ size 64
         proc (centered-spot size 610.0 2.0)
         imp (doto (ij.ImagePlus.) (.setProcessor proc))
         peaks (local-max-peaks imp)
         local-max (Util/copyPeakList peaks)
         data (imp-to-1d-double-array imp)
         width (.getWidth imp)
         height (.getHeight imp)
         multi (MultiFit. data peaks 1e-6 width height false)
         result (.getResults multi)
         foreground (.getForeground multi)
         residual (.getResidual multi)]
    (show-imp-overlay imp (peaks-update-overlay local-max result))
    (doseq [p local-max]
      (println (.getXWidth p) (.getYWidth p)))
    (visualize-double-buffer residual width height "Residual")
    (visualize-double-buffer foreground width height "Model")))

(defn test-single-centered-iter []
  (let [ size 64
         proc (centered-spot size 610.0 2.0)
         imp (doto (ij.ImagePlus.) (.setProcessor proc))
         peaks (local-max-peaks imp)
         local-max (Util/copyPeakList peaks)
         data (imp-to-1d-double-array imp)
         width (.getWidth imp)
         height (.getHeight imp)
         [result-peaks model residual] (do-fit imp peaks 1e-6 500 false)]
    (show-imp-overlay imp (peaks-update-overlay local-max result-peaks))
    (doseq [p local-max]
      (println (.getXWidth p) (.getYWidth p)))
    (visualize-double-buffer residual width height "Residual")
    (visualize-double-buffer model width height "Model")))


(test-fitting)

(comment
(defn fit-data
  [imp params &
   {:keys [background margin neighborhood new-peak-radius]
    :or {background -1 margin 10 neighborhood 5.0 new-peak-radius 1.0}}]
  (let [pad-size 10
        image (imp-to-1d-double-array imp)
        background (get :background params (imp-mean imp))
        cur-threshold (let [threshold (:threshold params)
                            iterations (:iterations params)]
                        (if (> iterations 4)
                          (* 4.0 threshold)
                          (double (* iterations threshold))))
        cutoff (+ background cur-threshold)]
    {:image image
     :background (get :background params (imp-mean imp))
     :cur-threshold cur-threshold
     :cutoff cutoff
     :find-max-radius 5
     :margin margin
     :neighborhood (* (:params sigma) neighborhood)
     :new-peak-radius new-peak-radius
     :pad-size (double pad-size)
     :peaks nil
     :proximity 5.0
     :residual image
     :sigma (:params sigma)
     :taken (int-array (alength image))
     :threshold (:threshold params)}))

  )

(defn -main [& args]
  (test-single-centered)
  (comment
  (let [ij (ij.ImageJ.)
        imp (.open (ij.io.Opener.) test-tiff)
        imp2 (.. imp getProcessor convertToFloat)]
    (do (.show imp)
        (.blurGaussian imp2 8)
        (.show (ij.ImagePlus. imp2))))))

