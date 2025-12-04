#' LI-6800 のデータを新しい葉面積で再計算する
#'
#' @param df racir::read_6800() で読み込んだデータフレーム
#' @param matched_area 新しい葉面積
#' @return 面積補正後のデータフレーム
#' @importFrom dplyr select mutate all_of any_of
#' @export
fixarea_6800 <- function(df, matched_area) {

  df |>
    dplyr::select(
      obs, A, Ci, CO2_r, CO2_s, Tleaf, E, gtc, Ca, S,
      H2O_r, H2O_s, CorrFact, Flow, Qin
    ) |>
    dplyr::mutate(
      matched_area = matched_area,
      A = Flow * CorrFact * (
        CO2_r - CO2_s * ((1000 - CorrFact * H2O_r) /
                           (1000 - CorrFact * H2O_s))
      ) / (100 * 10000 * matched_area),
      E = 1000 * CorrFact * Flow * (H2O_s - H2O_r) /
        (100 * 10000 * matched_area * (1000 - CorrFact / H2O_s)),
      Ci = -((gtc - E / 2) * Ca - A) / (gtc + E / 2),
      PPFD = Qin,
      gsc = A / (Ca - Ci),
      gsw = 1.6 * gsc,
      hhmmss = hhmmss
    )
}
