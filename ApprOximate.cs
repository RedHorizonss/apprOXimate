using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class ApprOximate : MonoBehaviour
{
  public TMP_InputField inputField;
  public TextMeshProUGUI outputText;
  public Toggle showFinalOnlyToggle;
  public CSVReader csvReader;

  public void AnalyzeFormula()
  {
    string formulaInput = inputField.text;
    var formula = ParseFormula(formulaInput);

    if (formula == null || formula.Count == 0)
    {
      outputText.text = "Invalid formula input.";
      return;
    }

    bool showFinalOnly = showFinalOnlyToggle != null && showFinalOnlyToggle.isOn;
    string output = PerformFullAnalysis(formula, includeIntermediate: !showFinalOnly, exportFinalComposition: false);
    outputText.text = output;
  }

  public void AnalyzeCSV()
  {
    if (csvReader == null)
    {
      Debug.LogError("CSVReader reference is missing.");
      return;
    }

    var formulas = csvReader.ReadCSV();
    if (formulas == null || formulas.Count == 0)
    {
      Debug.LogError("CSV file is empty or could not be read.");
      return;
    }

    bool showFinalOnly = showFinalOnlyToggle != null && showFinalOnlyToggle.isOn;

    List<string> results = new List<string>();
    foreach (var formulaInput in formulas)
    {
      var formula = ParseFormula(formulaInput);
      if (formula == null || formula.Count == 0)
      {
        results.Add($"Invalid formula: {formulaInput}");
      }
      else
      {
        string output = PerformFullAnalysis(formula, includeIntermediate: !showFinalOnly, exportFinalComposition: false);
        results.Add(output);
      }
    }

    outputText.text = string.Join("\n\n", results);
  }

  private string PerformFullAnalysis(Dictionary<string, double> formula, bool includeIntermediate, bool exportFinalComposition = false)
  {
    var localList1 = new Dictionary<string, List<(int Ox, int Red, double SRP)>>();
    foreach (var kvp in ChemicalData.List1)
      localList1[kvp.Key] = new List<(int Ox, int Red, double SRP)>(kvp.Value);

    double totalCharge;
    double previousChargeWithDelta = double.NaN;
    string output = includeIntermediate && !exportFinalComposition ? "" : "";

    string elementWithFinalLowestSRP = null;
    double finalDeltaOxMax = 0;
    double finalChargeWithDelta = 0;
    double elementQuantity = 0;
    int finalOx = 0, finalRed = 0;

    int iterationCount = 0;
    int maxIterations = 100;

    do
    {
      iterationCount++;
      if (iterationCount > maxIterations)
      {
        if (includeIntermediate && !exportFinalComposition)
          output += "<color=red>Maximum iterations reached. Stopping analysis to prevent freezing.</color>\n";
        break;
      }

      double sumNonM = 0;
      double sumM = 0;

      if (includeIntermediate && !exportFinalComposition)
        output += "Element Details:\n";

      foreach (var (element, quantity) in formula)
      {
        if (ChemicalData.List2[element] != -1)
        {
          int red = ChemicalData.List2[element];
          if (includeIntermediate && !exportFinalComposition)
            output += $"Element: {element}, Quantity: {Math.Round(quantity, 3)}, Red: {red}\n";
          sumNonM += quantity * red;
        }
        else
        {
          var data = localList1[element].FirstOrDefault(e => e.Red > 0);
          if (data != default)
          {
            if (includeIntermediate && !exportFinalComposition)
              output += $"Element: {element}, Quantity: {Math.Round(quantity, 3)}, Ox: {data.Ox}, Red: {data.Red}, SRP: {Math.Round(data.SRP, 3)}\n";
            sumM += quantity * data.Red;
          }
        }
      }

      totalCharge = Math.Round(sumNonM + sumM, 3);
      if (includeIntermediate && !exportFinalComposition)
      {
        output += $"Sum of (quant * Red) for non-M elements: {Math.Round(sumNonM, 3)}\n";
        output += $"Sum of (quant * Red) for M elements: {Math.Round(sumM, 3)}\n";
        output += $"Total Charge: {Math.Round(totalCharge, 3)}\n";
      }

      var lowestSrpElement = formula
          .Where(e => ChemicalData.List2[e.Key] == -1)
          .Select(e => new
          {
            Element = e.Key,
            Quantity = e.Value,
            Data = localList1[e.Key].FirstOrDefault(d => d.Red > 0)
          })
          .OrderBy(e => e.Data.SRP)
          .FirstOrDefault();

      if (lowestSrpElement != null)
      {
        double deltaOxMax = Math.Round((lowestSrpElement.Data.Ox - lowestSrpElement.Data.Red) * lowestSrpElement.Quantity, 3);
        double chargeWithDelta = Math.Round(deltaOxMax + totalCharge, 3);

        if (includeIntermediate && !exportFinalComposition)
        {
          output += $"M-element with the lowest SRP: {lowestSrpElement.Element}, SRP: {Math.Round(lowestSrpElement.Data.SRP, 3)}\n";
          output += $"{lowestSrpElement.Element} ΔOxMax = {Math.Round(deltaOxMax, 3)}\n";
          output += $"ΔOxMax + Total Charge = {Math.Round(chargeWithDelta, 3)}\n";
        }

        if (chargeWithDelta == previousChargeWithDelta)
        {
          if (includeIntermediate && !exportFinalComposition)
            output += "Stopping iteration as ΔOxMax + Total Charge did not change from the previous iteration.\n";
          break;
        }

        previousChargeWithDelta = chargeWithDelta;

        if (chargeWithDelta < 0)
        {
          totalCharge += deltaOxMax;
          var updatedData = localList1[lowestSrpElement.Element].FirstOrDefault(d => d.Red == lowestSrpElement.Data.Red);
          if (updatedData != default)
          {
            var newData = localList1[lowestSrpElement.Element].FirstOrDefault(d => d.Red == updatedData.Ox);
            if (newData == default)
              newData = (updatedData.Ox, updatedData.Ox, double.PositiveInfinity);

            localList1[lowestSrpElement.Element] = localList1[lowestSrpElement.Element].Where(d => d.Red != updatedData.Red).ToList();
            localList1[lowestSrpElement.Element].Add(newData);
          }
        }
        else
        {
          elementWithFinalLowestSRP = lowestSrpElement.Element;
          finalDeltaOxMax = deltaOxMax;
          finalChargeWithDelta = chargeWithDelta;
          elementQuantity = lowestSrpElement.Quantity;
          finalOx = lowestSrpElement.Data.Ox;
          finalRed = lowestSrpElement.Data.Red;
          break;
        }
      }
      else
      {
        if (includeIntermediate && !exportFinalComposition)
          output += "No valid M-element found for further iterations.\n";
        break;
      }

      if (includeIntermediate && !exportFinalComposition)
        output += "\n-----------------------------------\n\n";

    } while (true);

    double finalCharge = 0;
    foreach (var (element, quantity) in formula)
    {
      if (ChemicalData.List2[element] != -1)
      {
        finalCharge += quantity * ChemicalData.List2[element];
      }
      else if (element != elementWithFinalLowestSRP)
      {
        var finalData = localList1[element].FirstOrDefault(d => d.Red > 0);
        if (finalData != default)
          finalCharge += quantity * finalData.Red;
      }
    }

    if (finalChargeWithDelta == 0 || finalChargeWithDelta > 0)
    {
      double adjustedQuant = Math.Round(elementQuantity * (finalChargeWithDelta / finalDeltaOxMax), 3);
      double remainingQuant = Math.Round(elementQuantity - adjustedQuant, 3);

      if (adjustedQuant > 0)
        finalCharge += adjustedQuant * finalRed;

      if (remainingQuant > 0)
        finalCharge += remainingQuant * finalOx;
    }

    finalCharge = Math.Round(finalCharge, 3);
    string finalChargeSummary = $"FinalChargeBalance:{finalCharge}";

    if (exportFinalComposition)
    {
      List<string> segments = new List<string>();

      foreach (var (element, quantity) in formula)
      {
        if (ChemicalData.List2[element] != -1)
          segments.Add($"{element}:{ChemicalData.List2[element]}:{Math.Round(quantity, 3)}");
        else if (element != elementWithFinalLowestSRP)
        {
          var finalData = localList1[element].FirstOrDefault(d => d.Red > 0);
          if (finalData != default)
            segments.Add($"{element}:{finalData.Red}:{Math.Round(quantity, 3)}");
        }
      }

      if (finalChargeWithDelta == 0 || finalChargeWithDelta > 0)
      {
        double adjustedQuant = Math.Round(elementQuantity * (finalChargeWithDelta / finalDeltaOxMax), 3);
        double remainingQuant = Math.Round(elementQuantity - adjustedQuant, 3);

        if (adjustedQuant > 0)
          segments.Add($"{elementWithFinalLowestSRP}:{finalRed}:{adjustedQuant}");

        if (remainingQuant > 0)
          segments.Add($"{elementWithFinalLowestSRP}:{finalOx}:{remainingQuant}");
      }

      segments.Add(finalChargeSummary);
      return string.Join(";", segments);
    }
    else
    {
      foreach (var (element, quantity) in formula)
      {
        if (ChemicalData.List2[element] != -1)
          output += $"{Math.Round(quantity, 3)} {element} {ChemicalData.List2[element]}\n";
        else if (element != elementWithFinalLowestSRP)
        {
          var finalData = localList1[element].FirstOrDefault(d => d.Red > 0);
          if (finalData != default)
            output += $"{Math.Round(quantity, 3)} {element} {finalData.Red}\n";
        }
      }

      if (finalChargeWithDelta == 0 || finalChargeWithDelta > 0)
      {
        double adjustedQuant = Math.Round(elementQuantity * (finalChargeWithDelta / finalDeltaOxMax), 3);
        double remainingQuant = Math.Round(elementQuantity - adjustedQuant, 3);

        if (adjustedQuant > 0)
          output += $"{adjustedQuant} {elementWithFinalLowestSRP} {finalRed}\n";

        if (remainingQuant > 0)
          output += $"{remainingQuant} {elementWithFinalLowestSRP} {finalOx}\n";
      }

      string color = "green";
      if (Math.Abs(finalCharge) > 0.01) color = "red";
      else if (Math.Abs(finalCharge) > 0) color = "orange";

      output += $"<color={color}>{finalChargeSummary}</color>\n";
      return output;
    }
  }

  public void ExportFinalCompositionsToCSV()
  {
    if (csvReader == null)
    {
      Debug.LogError("CSVReader reference is missing.");
      return;
    }

    var formulas = csvReader.ReadCSV();
    if (formulas == null || formulas.Count == 0)
    {
      Debug.LogError("CSV file is empty or could not be read.");
      return;
    }

    List<string> results = new List<string>();
    foreach (var formulaInput in formulas)
    {
      var formula = ParseFormula(formulaInput);
      if (formula == null || formula.Count == 0)
        results.Add($"Invalid formula: {formulaInput}");
      else
        results.Add(PerformFullAnalysis(formula, includeIntermediate: false, exportFinalComposition: true));
    }

    string fileName = "final_compositions.csv";
    string filePath = Path.Combine(Application.persistentDataPath, fileName);

    try
    {
      File.WriteAllText(filePath, string.Join("\n", results));
      Debug.Log("Final compositions exported to: " + filePath);
    }
    catch (Exception e)
    {
      Debug.LogError("Error writing CSV file: " + e.Message);
    }
  }

  private Dictionary<string, double> ParseFormula(string formulaInput)
  {
    var regex = new Regex(@"([A-Z][a-z]*)(\d*\.?\d*\/?\d*)");
    var matches = regex.Matches(formulaInput);

    var formula = new Dictionary<string, double>();

    foreach (Match match in matches)
    {
      string element = match.Groups[1].Value;
      string quantityStr = match.Groups[2].Value;

      if (!ChemicalData.List2.ContainsKey(element))
        return null;

      double quantity = 1.0;

      if (!string.IsNullOrEmpty(quantityStr))
      {
        if (quantityStr.Contains("/"))
        {
          var fractionParts = quantityStr.Split('/');
          double numerator = double.Parse(fractionParts[0]);
          double denominator = double.Parse(fractionParts[1]);
          quantity = numerator / denominator;
        }
        else
        {
          quantity = double.Parse(quantityStr);
        }
      }

      if (quantity < 0.001)
        quantity = Math.Round(quantity, 6);

      quantity = Math.Round(quantity, 3);
      formula[element] = formula.ContainsKey(element) ? formula[element] + quantity : quantity;
    }

    return formula;
  }

  public void WriteOutputToCSV()
  {
    string fileName = "output.csv";
    string filePath = Path.Combine(Application.persistentDataPath, fileName);

    try
    {
      File.WriteAllText(filePath, outputText.text);
      Debug.Log("Output successfully written to: " + filePath);
    }
    catch (Exception e)
    {
      Debug.LogError("Error writing CSV file: " + e.Message);
    }
  }
}
