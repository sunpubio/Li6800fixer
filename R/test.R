#' LI-6800 のデータを新しい葉面積で再計算する
#'
#' @param df racir::read_6800() で読み込んだデータフレーム
#' @param matched_area 新しい葉面積
#' @return 面積補正後のデータフレーム
#' @importFrom dplyr select mutate all_of any_of
#' @export
fixarea_6800 <- function(df, matched_area) {

  df <- df[, !duplicated(names(df)), drop = FALSE]

  df |>
    dplyr::select(
      obs, A, Ci, CO2_r, CO2_s, Tleaf, E, gtc, Ca, S,
      H2O_r, H2O_s, CorrFact, Flow, Qin, TleafCnd, K, Pa, ΔPcham, gbw, hhmmss
    ) |>
    dplyr::mutate(
      matched_area = matched_area,
      A = Flow * CorrFact * (
        CO2_r - CO2_s * ((1000 - CorrFact * H2O_r) /
                           (1000 - CorrFact * H2O_s))
      ) / (100 * 10000 * matched_area),

      # ここは従来どおり（E は mmol ベースでOK）
      E = 1000 * CorrFact * Flow * (H2O_s - H2O_r) /
        (100 * 10000 * matched_area * (1000 - CorrFact / H2O_s)),

      Ci = -((gtc - E / 2) * Ca - A) / (gtc + E / 2),
      PPFD = Qin,

      # ★ gtw は E を mol に直して使う（E_mol = E / 1000）
      E_mol = E / 1000,

      gtw = E_mol * (
        1000 - ((1000 * 0.61365 * exp(17.502 * TleafCnd / (240.97 + TleafCnd)) /
                   (Pa + ΔPcham)) + H2O_s) / 2
      ) / (
        1000 * 0.61365 * exp(17.502 * TleafCnd / (240.97 + TleafCnd)) /
          (Pa + ΔPcham) - H2O_s
      ),

      gsw = 2 / (
        (1 / gtw - 1 / gbw) +
          sign(gtw) * sqrt(
            (1 / gtw - 1 / gbw)^2 +
              4 * K / ((K + 1)^2 * (2 * 1 / gtw * 1 / gbw - (1 / gbw)^2))
          )
      ),

      gtc = 1 / (
        (((K + 1) / (gsw / 1.6)) + 1 / (gbw / 1.37)) +
          (K / ((K + 1 / (gsw / 1.6)) + (K / (gbw / 1.37))))
      ),

      hhmmss = hhmmss
    )
}
